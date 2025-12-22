==========
Rattleback
==========

.. note::

    You can download this example as a Python script:
    :jupyter-download:script:`rattleback` or Jupyter notebook:
    :jupyter-download:notebook:`rattleback`.

Objectives
==========

- Show how to use ``symjit.compile_func`` as a (often) faster alternative to
  ``sympy.lambdify`` for compiling symbolic expressions.
- ``generate_ode_function`` from PyDy needs C - contiguous arrays as input,
  if generator='cython' is used. Depending on the method applied in
  ``scipy.integrate.solve_ivp``, it does not return C - contiguous arrays.
  Show how to handle this.
- Show a use of PyDy's visualization module to visualize the results of
  a simulation.

Description
===========

This simulation is taken from this article:
https://www.sciencedirect.com/science/article/abs/pii/0020746282900178
The equations of motion are set up using Kane's method from SymPy's
``physics.mechanics`` module. The ingenious way to get the contact points of
the rattleback with the plane is from this paper.

Notes:
======

- The numerical integration seems very sensitive:

  - method=DOP853, rtol=1.e-5, atol=1.e-6 seems the minimum to get a
    qualitatively reasonable result.
  - changing C = :math:`1.7 \cdot 10^{-3}` (one of the moments of inertia)
    to C = :math:`1.65 \cdot 10^{-3}`
    or C = :math:`1.75 \cdot 10^{-3}` makes the integration fail.
  - changing the initial angular velocity :math:`u_3` from -1.0 to -1.25
    makes the integration fail. The paper gives -5.0.
  - An indication of the difficulty is that to integrate an interval of 5 sec,
    ``solve_ivp`` needs more than 800,000 function evaluations.

- The total energy does not decrease as expected, maybe due to numerical
  problems



.. jupyter-execute::

    import numpy as np
    import sympy as sm
    import symjit
    import sympy.physics.mechanics as me
    import matplotlib.pyplot as plt
    from copy import deepcopy
    from scipy.integrate import solve_ivp
    from pydy.codegen.ode_function_generators import generate_ode_function
    from pydy.viz.shapes import Cylinder, Sphere, Box
    from pydy.viz.scene import Scene
    from pydy.viz.visualization_frame import VisualizationFrame


Set up the Equations of Motion
------------------------------

.. jupyter-execute::

    N, R = sm.symbols('N, R', cls=me.ReferenceFrame)
    O = me.Point('O')
    Oe = me.Point('Oe')
    t = me.dynamicsymbols._t
    O.set_vel(N, 0)

``Ro`` is the mass center of the rattleback, ``S`` is the contact point with
the plane.

.. jupyter-execute::

    S, Ro = sm.symbols('S, Ro', cls=me.Point)

Rotation angles and speeds of ``R`` w.r.t. ``N``

.. jupyter-execute::

    q1, q2, q3, u1, u2, u3 = me.dynamicsymbols('q1, q2, q3, u1, u2, u3')
    x1, x2, x3, ux1, ux2, ux3 = me.dynamicsymbols('x1, x2, x3, ux1, ux2, ux3')
    xe, ye, ze, uxe, uye, uze = me.dynamicsymbols('xe, ye, ze, uxe, uye, uze')

symjit cannot use dynamicsymbols, so sympy symbols for replacement in the
compilation are defined.

.. jupyter-execute::

    qs1, qs2, qs3 = sm.symbols('qs1, qs2, qs3')
    us1, us2, us3 = sm.symbols('us1, us2, us3')
    xs1, xs2, xs3 = sm.symbols('xs1, xs2, xs3')
    uxs1, uxs2, uxs3 = sm.symbols('uxs1, uxs2, uxs3')

Some physical parameters.

.. jupyter-execute::

    h, g, m, A, B, C, D = sm.symbols('h, g, m, A, B, C, D')
    a, b, c = sm.symbols('a, b, c')
    friktion = sm.symbols('friktion')

Set up the geometry and the kinematics.

.. jupyter-execute::

    R.orient_body_fixed(N, (q1, q2, q3), 'XYZ')
    rot = R.ang_vel_in(N)
    R.set_ang_vel(N, u1 * R.x + u2 * R.y + u3 * R.z)
    rot1 = R.ang_vel_in(N)

``Oe`` is the geometric center of the rattleback. The dynamics of the system
does not depend on it.

.. jupyter-execute::

    Oe.set_pos(O, xe * N.x + ye * N.y + ze * N.z)
    Oe.set_vel(N, uxe * N.x + uye * N.y + uze * N.z)

``Ro`` is the mass center of the rattleback.

.. jupyter-execute::

    Ro.set_pos(Oe, -h * R.z)

The coordinates and the speeds of the contact point ``S`` in ``R`` will be
calculated numerically during the integration.

.. jupyter-execute::

    S.set_pos(Oe, x1*R.x + x2*R.y + x3*R.z)
    S.set_vel(R, ux1*R.x + ux2*R.y + ux3*R.z)

No slip condition: ``S`` is at rest in ``N`` momentarily.

.. jupyter-execute::

    Ro.set_vel(N, R.ang_vel_in(N).cross(Ro.pos_from(S)))

Define the rigid body of the rattleback.

.. jupyter-execute::

    inert = me.inertia(R, A, B, C, D, 0, 0)
    rattleback = me.RigidBody('Rattleback', Ro, R, m, (inert, Ro))

Finish setting up Kane's equations.

.. jupyter-execute::

    bodies = [rattleback]

    forces = [
        (Ro, -m * g * N.z),
        (R, -friktion * R.ang_vel_in(N))
    ]

    q_ind = [q1, q2, q3]
    u_ind = [u1, u2, u3]

    kd = sm.Matrix([(rot - rot1).dot(uv) for uv in N])

    kane = me.KanesMethod(N, q_ind, u_ind, kd_eqs=kd)
    fr, frstar = kane.kanes_equations(bodies, forces)
    mass_matrix = kane.mass_matrix
    force = me.msubs(kane.forcing, {x1.diff(t): ux1,
                                    x2.diff(t): ux2,
                                    x3.diff(t): ux3})

Print some information about the equations of motion.

.. jupyter-execute::

    print('Mass matrix dynamic symbols:', me.find_dynamicsymbols(mass_matrix))
    print('Mass matrix free symbols:', mass_matrix.free_symbols, '\n')
    print('Force dynamic symbols:', me.find_dynamicsymbols(force))
    print('Force free symbols:', force.free_symbols)
    print(f"Force contains {sm.count_ops(force)} operations.")


Generate and Compile the Equations of Motion
--------------------------------------------

.. jupyter-execute::

    pL = [m, g, a, b, c, h, A, B, C, D, friktion, xs1, xs2, xs3, uxs1, uxs2,
          uxs3]

    symbol_dict_1 = {
        x1: xs1,
        x2: xs2,
        x3: xs3,
        ux1: uxs1,
        ux2: uxs2,
        ux3: uxs3,
    }

    symbol_dict_2 = {
        q1: qs1,
        q2: qs2,
        q3: qs3,
        u1: us1,
        u2: us2,
        u3: us3,
    }

    specified = None
    constants = np.array(pL)
    loesung = sm.solve(kd, [q_ind[i].diff(t) for i in range(3)])

The solution for the kinematic equations must be sorted so that it corresponds
to kane.q

.. jupyter-execute::

    schluessel = [i.diff(t) for i in q_ind]
    kin_eqs_solved = sm.Matrix([loesung[i] for i in schluessel])

Create the rhs function.

.. jupyter-execute::

    force = me.msubs(force, symbol_dict_1)
    mass_matrix = me.msubs(mass_matrix, symbol_dict_1)

    rhs_gen = generate_ode_function(
        force,
        kane.q,
        kane.u,
        constants=constants,
        mass_matrix=mass_matrix,
        specifieds=specified,
        coordinate_derivatives=kin_eqs_solved,  # rhs of kin. diff. equations
        generator='cython',
        linear_sys_solver='sympy',
        constants_arg_type='array',
        specifieds_arg_type='array',
    )


Get the position of the contact point, using Kane's explicit solution given in
the paper mentioned above. As they are evaluated in every call to the r.h.s.
during numerical integration, compile with ``symjit`` for speed.

.. jupyter-execute::

    pL_pos = [m, g, a, b, c, h, A, B, C, D, friktion]

    epsilon = -sm.sqrt((a * N.z.dot(R.x))**2 + (b * N.z.dot(R.y))**2 +
                       (c * N.z.dot(R.z))** 2)

    x11 = (a**2 * N.z.dot(R.x) / epsilon)
    x21 = (b**2 * N.z.dot(R.y) / epsilon)
    x31 = (c**2 * N.z.dot(R.z) / epsilon)

    symbol_dict = symbol_dict_1 | symbol_dict_2

    x1_sub = x11.subs(symbol_dict)
    x2_sub = x21.subs(symbol_dict)
    x3_sub = x31.subs(symbol_dict)

    pos_kane = [x1_sub, x2_sub, x3_sub]
    pos_jit = symjit.compile_func([qs1, qs2, qs3, us1, us2, us3], pos_kane,
                              params=pL_pos)


Get the speed of (subsequent) contact points, compile with ``symjit``.

.. jupyter-execute::

    pL_speed = [m, g, a, b, c, h, A, B, C, D, friktion, xs1, xs2, xs3]

    pos_S = x11 * R.x + x21 * R.y + x31 * R.z
    pos_S_dt = pos_S.diff(t, R)
    coord_dict = {q1.diff(t): u1, q2.diff(t): u2, q3.diff(t): u3}

    x11dt = pos_S_dt.dot(R.x).subs(coord_dict)
    x21dt = pos_S_dt.dot(R.y).subs(coord_dict)
    x31dt = pos_S_dt.dot(R.z).subs(coord_dict)

    x11dt_sub = x11dt.subs(symbol_dict)
    x21dt_sub = x21dt.subs(symbol_dict)
    x31dt_sub = x31dt.subs(symbol_dict)

    speed_kane = [x11dt_sub, x21dt_sub, x31dt_sub]
    speed_jit = symjit.compile_func([qs1, qs2, qs3,
                                     us1, us2, us3],
                                    speed_kane,
                                    params=pL_speed)

While the speed of ``S`` has no ``N.z`` component, the position of ``S``
in ``N.z`` is not constant, as per the set up to get Kane's equations. To
compensate for this, ``ze_wert`` is defined.
Needed later for the potential energy.

.. jupyter-execute::

    ze_wert = sm.solve(S.pos_from(O).dot(N.z), ze)[0]


Show numerically that

- the speed of ``S`` has zero ``N.z`` component within numerical accuracy
- ``ze_wert`` 'compensates' ``S.pos_from(Oe)`` within numerical accuracy

.. jupyter-execute::

    kd_solved = sm.solve(kd, (q1.diff(t), q2.diff(t), q3.diff(t)))

    test = (x11dt * R.x + x21dt * R.y + x31dt * R.z).dot(N.z)
    test = test.subs(kd_solved)
    test_lam = sm.lambdify((q1, q2, q3, u1, u2, u3, a, b, c), test, cse=True)

    S_coord_N =[(x11 * R.x + x21 * R.y + x31 * R.z).dot(uv) for uv in N]
    S_coord_N_lam = sm.lambdify((q1, q2, q3, a, b, c), S_coord_N, cse=True)
    S_coords_R_lam = sm.lambdify((q1, q2, q3, a, b, c), [x11, x21, x31],
                                 cse=True)
    Sz_wert_lam = sm.lambdify((q1, q2, q3, x1, x2, x3, a, b, c), ze_wert,
                              cse=True)

    a1 = 1
    b1 = 2
    c1 = 3
    zaehler = 0
    for i in range(50):
        q1v = np.random.rand() * 2 * np.pi
        q2v = np.random.rand() * 2 * np.pi
        q3v = np.random.rand() * 2 * np.pi
        u1v = np.random.randn() * 10
        u2v = np.random.randn() * 10
        u3v = np.random.randn() * 10
        x11v, x21v, x31v = S_coords_R_lam(q1v, q2v, q3v, a1, b1, c1)
        Sz_wert1 = Sz_wert_lam(q1v, q2v, q3v, x11v, x21v, x31v, a1, b1, c1)
        res = test_lam(q1v, q2v, q3v, u1v, u2v, u3v, a1, b1, c1)
        Sz = S_coord_N_lam(q1v, q2v, q3v, a1, b1, c1)[2]
        if np.abs(res) > 1.e-13:
            zaehler += 1
            print(f"Result {i + 1}: {res}")
        if np.abs(Sz + Sz_wert1) > 1.e-13:
            zaehler += 1
            print(f"Sz check {i + 1}: {Sz} + {Sz_wert1} = {Sz}")

    if zaehler == 0:
        print("All tests passed successfully.")


Numerical Integration
---------------------

Seems quite difficult to integrate numerically. If rtol and atol are too
large, the integration may finish successfully, but the results will be
qualitatively different.


Input values, taken from the paper mentioned above.

.. jupyter-execute::

    m1 = 1.                          # mass of the ball
    g1 = 9.81                        # gravitational acceleration
    a1 = 0.2                         # semi axes of the ellipsoid
    b1 = 0.03                        # semi axes of the ellipsoid
    c1 = 0.02                        # semi axes of the ellipsoid
    h1 = 0.01                        # distance geometric center to mass center
    A1 = 2.e-4
    B1 = 1.6e-3
    C1 = 1.7e-3
    D1 = -2.0e-5
    friktion1 = 1.e-4
    q11 = np.deg2rad(0.5)        # initial orientation angles
    q21 = np.deg2rad(0.5)
    q31 = np.deg2rad(0.0)

    u11 = 0.0                     # initial angular velocities
    u21 = 0.0

In the paper :math:`u_{31}` is set to -5 rad/sec. Here, if it is decreased to
less than -1.0 the numerical integration will fail.

.. jupyter-execute::

    u31 = -1.0

.. jupyter-execute::

    xs11, xs21, xs31, uxs11, uxs21, uxs31 = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    pL_vals = [m1, g1, a1, b1, c1, h1, A1, B1, C1, D1, friktion1,
               xs11, xs21, xs31, uxs11, uxs21, uxs31]

    y0 = [q11, q21, q31,
          u11, u21, u31]

    interval = 5.0
    schritte = int(1000 * interval)
    times = np.linspace(0., interval, schritte)

    t_span = (0., interval)
    zaehler = 0
    args_list = []
    time_list = []

Define the right hand side function for the integrator. In each step the
position and speed of the contact point S are calculated and stored in
args_list.

.. jupyter-execute::

    def gradient(t, y, args):
        global zaehler
        zaehler += 1

        # Calculate the location of the contact point S
        args1 = [args[i] for i in range(len(args)-6)]
        args[-6], args[-5], args[-4] = pos_jit(*y, *args1)

        # Calculate the speed of the contact point S
        args2 = [args[i] for i in range(len(args)-3)]
        args[-3], args[-2], args[-1] = speed_jit(*y, *args2)

        # Store args for later plotting, energy, etc.
        time_list.append(t)
        args_list.append(deepcopy(args))
        args = np.array(args)
        # Ensure y is C - contiguous
        y = np.ascontiguousarray(y)
        rhs = rhs_gen(y, t, args)
        return rhs

Find the solution.

.. jupyter-execute::

    resultat1 = solve_ivp(gradient, t_span, y0, t_eval=times, args=(pL_vals,),
                          method='DOP853',
                          rtol=1.e-12,
                          atol=1.e-13)

    resultat = resultat1.y.T

    print(resultat1.message)
    print('Shape of resultat', resultat.shape)
    print("To numerically integrate an intervall of {:.3f} sec the routine "
          "cycled {:,} times".format(interval, resultat1.nfev))


Plot the Rotation Angles and Speeds
-----------------------------------

.. jupyter-execute::

    fig, ax = plt.subplots(2, 1, figsize=(8, 3), constrained_layout=True,
                           sharex=True)
    for i in range(3):
        ax[0].plot(resultat1.t, np.rad2deg(resultat[:, i]), label=f'q{i+1}')
        ax[1].plot(resultat1.t, resultat[:, i+3], label=f'u{i+1}')
        ax[1].axhline(0, color='k', lw=0.5, ls='--')

    ax[0].set_ylabel('[deg]')
    ax[0].set_title('Generalized Coordinates')
    ax[1].set_ylabel('[rad/s]')
    ax[1].set_title('Angular Velocities')
    ax[1].set_xlabel('Time [s]')
    ax[0].legend()
    _ = ax[1].legend()

Plot the Energy
---------------

It does not behave as expected, it should drop due to friction being present.
Maybe due to numerical inaccuracies.

``args_list`` and ``time_list`` have many more entries than ``resultat``
(850,000 vs. 5000) so the entries in time_list (and hence in args_list) that
correspond closest to the times in ``times`` must be found. Errors here may
also contribute to the energy not being perfectly conserved.

.. jupyter-execute::

    time_list = np.array(time_list)
    idxs = np.searchsorted(time_list, times)
    idxs_clipped = np.clip(idxs, 1, len(time_list) - 1)  # avoid out-of-bounds

    left = idxs_clipped - 1
    right = idxs_clipped

    idxs_final = np.where(
        np.abs(times - time_list[left]) <= np.abs(times - time_list[right]),
        left, right)

Define and compile the energy expressions.

.. jupyter-execute::

    qL = [q1, q2, q3, u1, u2, u3]

    pot_energy = m * g * (Ro.pos_from(Oe).dot(N.z).subs(symbol_dict_1) +
                          ze_wert.subs(symbol_dict_1))
    kin_energy = rattleback.kinetic_energy(N).subs(symbol_dict_1)
    pot_lam = sm.lambdify(qL + pL, pot_energy, cse=True)
    kin_lam = sm.lambdify(qL + pL, kin_energy, cse=True)

    pot_np = np.empty(resultat.shape[0])
    kin_np = np.empty(resultat.shape[0])
    total_np = np.empty(resultat.shape[0])

    for i in range(resultat.shape[0]):
        pL1_vals = args_list[idxs_final[i]]
        pot_np[i] = pot_lam(*resultat[i], *pL1_vals)
        kin_np[i] = kin_lam(*resultat[i], *pL1_vals)
        total_np[i] = pot_np[i] + kin_np[i]

    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 2), constrained_layout=True)
    ax2.plot(resultat1.t, pot_np, label='Potential Energy')
    ax2.plot(resultat1.t, kin_np, label='Kinetic Energy')
    ax2.plot(resultat1.t, total_np, label='Total Energy')
    ax2.set_ylabel('[Nm]')
    ax2.set_title('Energies')
    ax2.set_xlabel('Time [s]')
    _ = ax2.legend()

Coordinates and speeds of the Contact Point S as seen in R and N.

.. jupyter-execute::

    bezeichnung = ['x1', 'x2', 'x3', 'ux1', 'ux2', 'ux3']
    fig, ax = plt.subplots(4, 1, figsize=(8, 6), sharex=True,
                           constrained_layout=True)
    for i in range(3):
        ax[0].plot([time_list[idxs_final[j]] for j in range(resultat.shape[0])],
                   [args_list[idxs_final[j]][-6 + i]
                    for j in range(resultat.shape[0])], label=bezeichnung[i])
    ax[0].set_ylabel('[m]')
    ax[0].set_title('Position of Contact Point as seen in R')
    _ = ax[0].legend()

    for i in range(3):
        ax[2].plot([time_list[idxs_final[j]] for j in range(resultat.shape[0])],
                   [args_list[idxs_final[j]][-3 + i]
                    for j in range(resultat.shape[0])],
            label=bezeichnung[i + 3])
        ax[2].set_ylabel('[m/s]')
        ax[2].set_title('Speed of Contact Point as seen in R')
        _ = ax[2].legend()

    repl_dict = symbol_dict_1 | {q1.diff(t): u1,
                                 q2.diff(t): u2,
                                 q3.diff(t): u3,
                                 sm.Derivative(xs1, t): uxs1,
                                 sm.Derivative(xs2, t): uxs2,
                                 sm.Derivative(xs3, t): uxs3,
                                }

    S_N = S.pos_from(Oe)
    S_N = [S_N.dot(N.x).subs(repl_dict),
           S_N.dot(N.y).subs(repl_dict),
           S_N.dot(N.z).subs(repl_dict) + ze_wert.subs(repl_dict)]

    Sdt_N = (S.pos_from(Oe).diff(t, R))
    Sdt_N = [Sdt_N.dot(N.x).subs(repl_dict),
             Sdt_N.dot(N.y).subs(repl_dict),
             Sdt_N.dot(N.z).subs(repl_dict)]

    S_N_lam = sm.lambdify(qL + pL, S_N, cse=True)
    Sdt_N_lam = sm.lambdify(qL + pL, Sdt_N, cse=True)

    bezeichnung = ['S_N_x', 'S_N_y', 'S_N_z', 'Sdt_N_x', 'Sdt_N_y', 'Sdt_N_z']

    for i in range(3):
        ax[1].plot([time_list[idxs_final[j]] for j in range(resultat.shape[0])],
                   [S_N_lam(*resultat[j], *args_list[idxs_final[j]])
                    [i] for j in range(resultat.shape[0])],
                   label=bezeichnung[i])
    ax[1].set_ylabel('[m]')
    ax[1].set_title('Position of Contact Point as seen in N')
    _ = ax[1].legend()

    for i in range(3):
        ax[3].plot([time_list[idxs_final[j]] for j in range(resultat.shape[0])],
                    [Sdt_N_lam(*resultat[j], *args_list[idxs_final[j]])[i]
                    for j in range(resultat.shape[0])],
                   label=bezeichnung[i + 3])
    ax[3].set_ylabel('[m/s]')
    ax[3].set_title('Speed of Contact Point as seen in N')
    ax[3].set_xlabel('Time [s]')
    _ = ax[3].legend()


The ``N.z`` component of the speed of ``S`` should be zero within numerical
accuracy. The ``N.z`` component of the position of ``S`` should be compensated
by ``ze_wert`` to zero within numerical accuracy. While the position seems to
be very accurate, the speed is less so. This may affect the kinetic energy.

.. jupyter-execute::

    Nz_deviation = [S_N_lam(*resultat[j], *args_list[idxs_final[j]])
                            [2] for j in range(resultat.shape[0])]
    Nzdt_deviation = [Sdt_N_lam(*resultat[j], *args_list[idxs_final[j]])
                                [2] for j in range(resultat.shape[0])]
    print(f"Maximal N.z deviation of S: {np.max(np.abs(Nz_deviation)):.2e}",
          "ideally should be zero")
    print(f"Maximal d(N.z)/dt deviation of S: "
          f"{np.max(np.abs(Nzdt_deviation)):.2e}"
          " ideally should be zero")

Visualization
-------------

.. jupyter-execute::

    point1, point2, point3 = sm.symbols('point1, point2, point3', cls=me.Point)
    point1.set_pos(Oe, a/2 * R.x + c/2 * R.z)
    point2.set_pos(Oe, a * R.x)
    farben = ['red', 'green', 'blue']

    viz_frames = []

``groesse`` is an empirical value to make the animation a good size.

.. jupyter-execute::

    groesse = 2

As the axes of pythreejs are fixed, and the equations of motion were set up
differently, they are rotated to match pythreejs' axes.

.. jupyter-execute::

    B = me.ReferenceFrame('B')
    B.orient_body_fixed(N, (sm.pi/2, 0, 0), 'XYZ')

Start the animation.

.. jupyter-execute::

    rattle_shape = Box(name='rattle',
                       width=a * groesse,
                       height=b * groesse,
                       depth=c * groesse,
                       color='grey')

    viz_frames.append(VisualizationFrame('rattle_frame',
                                         rattleback,
                                         rattle_shape))

    for i, point in enumerate([point1, point2]):
        point_shape = Sphere(name='point{}'.format(i), radius=0.01 * groesse,
                             color=farben[i])
        viz_frames.append(VisualizationFrame('point_frame{}'.format(i),
                                             R,
                                             point,
                                             point_shape))

    scene = Scene(B, Oe, *viz_frames)

    scene.times = times
    pL_vals_adjusted = [val * groesse for val in pL_vals]
    scene.constants = dict(zip(pL, pL_vals_adjusted))
    scene.states_symbols = q_ind + u_ind
    scene.states_trajectories = resultat
    scene.display_jupyter(axes_arrow_length=20)
