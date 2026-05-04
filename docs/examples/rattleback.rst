==========
Rattleback
==========

.. note::

    You can download this example as a Python script:
    :jupyter-download:script:`rattleback` or Jupyter notebook:
    :jupyter-download:notebook:`rattleback`.

Objectives
==========

- Show how to define a special ODE solver.
- Show a use of PyDy's visualization module to visualize the results of
  a simulation.

Description
===========

This simulation is taken from this article:
https://www.sciencedirect.com/science/article/abs/pii/0020746282900178
The equations of motion are set up using Kane's method from SymPy's
``physics.mechanics`` module. The ingenious way to get the contact points of
the rattleback with the plane is from the paper.

Notes:
======

- The numerical integration seems quite difficult. Even with ''method=DOP853'',
  which is a high order method, the standard values of ``atol`` and ``rtol``
  are not sufficient to get a qualitatively reasonable result.
- It seems to depend critically on the constants, also taken from the a.m.
  paper. (Maybe this is the reason rattleback are not seen often in nature)


.. jupyter-execute::

    import numpy as np
    import sympy as sm
    import sympy.physics.mechanics as me
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from scipy.integrate import solve_ivp
    from pydy.system import System
    from pydy.viz.shapes import Sphere, Box
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

Rotation angles and speeds.

.. jupyter-execute::

    q1, q2, q3, u1, u2, u3 = me.dynamicsymbols('q1, q2, q3, u1, u2, u3')
    x1, x2, x3, ux1, ux2, ux3 = me.dynamicsymbols('x1, x2, x3, ux1, ux2, ux3')
    xe, ye, ze, uxe, uye, uze = me.dynamicsymbols('xe, ye, ze, uxe, uye, uze')


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

``Oe`` is the geometric center of the rattleback.

.. jupyter-execute::

    Oe.set_pos(O, xe * N.x + ye * N.y + ze * N.z)

``Ro`` is the mass center of the rattleback.

.. jupyter-execute::

    Ro.set_pos(Oe, -h * R.z)

The coordinates and the speeds of the contact point ``S`` in ``R`` will be
dependent coordinates and speeds.


.. jupyter-execute::

    S.set_pos(Oe, x1*R.x + x2*R.y + x3*R.z)
    S.set_vel(N, ux1*R.x + ux2*R.y + ux3*R.z)

No slip condition: ``S`` is at rest in ``N`` momentarily is used to get the
velocity of Oe.

.. jupyter-execute::

    Oe.set_vel(N, R.ang_vel_in(N).cross(Oe.pos_from(S)))
    Ro.set_vel(N, R.ang_vel_in(N).cross(Ro.pos_from(Oe)) + Oe.vel(N))

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

Calculate the position of the contact point as per the a.m. paper.

.. jupyter-execute::

    epsilon = -sm.sqrt((a * N.z.dot(R.x))**2 + (b * N.z.dot(R.y))**2 +
                       (c * N.z.dot(R.z))** 2)

    config_constr = sm.Matrix([
            x1 - (a**2 * N.z.dot(R.x) / epsilon),
            x2 - (b**2 * N.z.dot(R.y) / epsilon),
            x3 - (c**2 * N.z.dot(R.z) / epsilon)
        ])

Get the speed constraints for Oe.

.. jupyter-execute::

    speed_constr_S = config_constr.diff(t)

    speed_constr_Oe = sm.Matrix([
            uxe - Oe.vel(N).dot(N.x),
            uye - Oe.vel(N).dot(N.y),
            uze - Oe.vel(N).dot(N.z),
    ])

Combine the speed constraints.

.. jupyter-execute::

    speed_constr = sm.Matrix([*speed_constr_S, *speed_constr_Oe])

Get Kane's equations.

.. jupyter-execute::

    q_ind = [q1, q2, q3, x1, x2, x3, xe, ye, ze]
    u_ind = [u1, u2, u3]
    u_dep = [ux1, ux2, ux3, uxe, uye, uze]

    kd = sm.Matrix([*[((rot - rot1).dot(uv)) for uv in N],
                    x1.diff(t) - ux1,
                    x2.diff(t) - ux2,
                    x3.diff(t) - ux3,
                    xe.diff(t) - uxe,
                    ye.diff(t) - uye,
                    ze.diff(t) - uze
                  ])

    kane = me.KanesMethod(
        N,
        q_ind,
        u_ind,
        u_dependent=u_dep,
        kd_eqs=kd,
        velocity_constraints=speed_constr,
    )
    fr, frstar = kane.kanes_equations(bodies, forces)


Set up an ODE solver for increased accuracy needed for this system.

.. jupyter-execute::

    def ode_solver(func, y0, times, args=(), **kwargs):
        res = solve_ivp(lambda t, y, *args: func(y, t, *args),
                       (times[0], times[-1]),
                        y0,
                        args=args,
                        t_eval=times,
                        method='DOP853',
                        atol=1.e-12,
                        rtol=1.e-12,
                        **kwargs)
        return res.y.T

Set up an instance of System.

.. jupyter-execute::

    sys = System(kane, ode_solver=ode_solver)

Define the constants of the system, taken from the paper mentioned above.

.. jupyter-execute::

    sys.constants = {
    m: 1.0,                         # mass of the ball
    g: 9.81,                        # gravitational acceleration
    a: 0.2,                         # semi axes of the ellipsoid
    b: 0.03,                        # semi axes of the ellipsoid
    c: 0.02,                        # semi axes of the ellipsoid
    h: 0.01,                        # distance geometric center to mass center
    A: 2.e-4,                       # inertia parameters A, B, C, D
    B: 1.6e-3,
    C: 1.7e-3,
    D: -2.0e-5,
    friktion: 1.e-4,                # friction coefficient
    }

    print('Constants of the system:')
    sys.constants


Set the initials conditions of the independent coordinates and speeds.
Values taken from the paper mentioned above.

.. jupyter-execute::

    sys.initial_conditions = {
    q1: np.deg2rad(0.5),        # initial orientation angles
    q2: np.deg2rad(0.5),
    q3: np.deg2rad(0.0),
    u1: 0.0,                     # initial angular velocities
    u2: 0.0,
    u3: -1.0,
    }


Get the position and speed of the contact point for consistent initial
conditions. The position of Oe ist set such that S is at O initially,
and the speed of Oe is set such that S is at rest in N initially.

.. jupyter-execute::

    pL_pos = [key for key in sys.constants.keys()]
    pL_vals = [sys.constants[key] for key in pL_pos]
    qL = [q1, q2, q3, u1, u2, u3]
    qL_vals = [sys.initial_conditions[key] for key in qL]

Pos of S in R as per D Levinson's paper

.. jupyter-execute::

    S_pos = sm.solve(config_constr, (x1, x2, x3))
    S_pos_R = [S_pos[x1], S_pos[x2], S_pos[x3]]
    S_pos_R_lam = sm.lambdify(qL + pL_pos, S_pos_R, cse=True)
    x1_1, x2_1, x3_1 = S_pos_R_lam(*qL_vals, *pL_vals)
    dict_S_pos = {x1: x1_1, x2: x2_1, x3: x3_1}

Speed of S in R.

.. jupyter-execute::

    speed_constr_S = speed_constr_S.subs({i.diff(t): j
                                      for i, j in zip(q_ind, u_ind + u_dep)})
    S_speed_R = sm.solve(speed_constr_S, (ux1, ux2, ux3))
    S_speed_R = [S_speed_R[ux1].subs(dict_S_pos),
                 S_speed_R[ux2].subs(dict_S_pos),
                 S_speed_R[ux3].subs(dict_S_pos)]
    S_speed_R_lam = sm.lambdify(qL + pL_pos, S_speed_R, cse=True)
    ux1_1, ux2_1, ux3_1 = S_speed_R_lam(*qL_vals, *pL_vals)
    dict_S_speed_R = {ux1: ux1_1, ux2: ux2_1, ux3: ux3_1}

Initial conditons for Oe in N, so that S is at O initially.

.. jupyter-execute::


    S_pos_N = sm.Matrix([S.pos_from(O).dot(N.x), S.pos_from(O).dot(N.y),
                         S.pos_from(O).dot(N.z)])

    Oe_pos_N = sm.solve(S_pos_N, (xe, ye, ze))
    Oe_pos_N = [Oe_pos_N[xe].subs(S_pos), Oe_pos_N[ye].subs(S_pos),
                Oe_pos_N[ze].subs(S_pos)]
    Oe_pos_N_lam = sm.lambdify(qL + pL_pos, Oe_pos_N, cse=True)
    xe1, ye1, ze1 = Oe_pos_N_lam(*qL_vals, *pL_vals)
    dict_Oe_pos_N = {xe: xe1, ye: ye1, ze: ze1}

Initial conditions for Oe speed in N, so that S is at rest initially.

.. jupyter-execute::

    Oe_speed_N = sm.solve(speed_constr_Oe, (uxe, uye, uze))
    Oe_speed_N = [me.msubs(Oe_speed_N[uxe], dict_Oe_pos_N, dict_S_pos),
                  me.msubs(Oe_speed_N[uye], dict_Oe_pos_N, dict_S_pos),
                  me.msubs(Oe_speed_N[uze], dict_Oe_pos_N, dict_S_pos)]
    Oe_speed_N_lam = sm.lambdify(qL + pL_pos, Oe_speed_N)
    uxe1, uye1, uze1 = Oe_speed_N_lam(*qL_vals, *pL_vals)
    dict_Oe_speed_N = {uxe: uxe1, uye: uye1, uze: uze1}

Print the initial conditions.

.. jupyter-execute::

    sys.initial_conditions = {**sys.initial_conditions, **dict_S_pos,
                          **dict_S_speed_R, **dict_Oe_pos_N,
                          **dict_Oe_speed_N}

    print('Initial conditions:')
    dict_clean = {k: float(v) for k, v in sys.initial_conditions.items()}
    dict_clean


Numerical Integration
---------------------

.. jupyter-execute::

    sys.generate_ode_function(generator='cython', linear_sys_solver='numpy')

    interval = 5.0
    schritte = int(1000 * interval)
    sys.times = np.linspace(0., interval, schritte)
    times = sys.times

    resultat = sys.integrate()
    print('Shape of resultat', resultat.shape)


Plot Some Results
-----------------

.. jupyter-execute::

    qL1 = q_ind + u_ind + u_dep
    bezeichnung = [str(i) for i in qL1]

    S_vel_N = S.vel(N).subs({i.diff(t): j for i, j in zip(q_ind,
                                                          u_ind + u_dep)})
    S_vel_N = sm.Matrix([S_vel_N.dot(N.x), S_vel_N.dot(N.y), S_vel_N.dot(N.z)])
    S_vel_N_lam = sm.lambdify(qL1 + pL_pos, S_vel_N, cse=True)
    S_vel_N_np = np.zeros((resultat.shape[0], 3))
    for i in range(resultat.shape[0]):
        S_vel_N_np[i] = S_vel_N_lam(*[resultat[i, j]
                                      for j in range(resultat.shape[1])],
                                    *pL_vals).squeeze()

    S_pos_N = [S.pos_from(O).dot(N.x), S.pos_from(O).dot(N.y),
               S.pos_from(O).dot(N.z)]
    S_pos_N_lam = sm.lambdify(qL1 + pL_pos, S_pos_N, cse=True)
    S_pos_np = np.zeros((resultat.shape[0], 3))
    for i in range(resultat.shape[0]):
        S_pos_np[i] = S_pos_N_lam(*[resultat[i, j]
                                    for j in range(resultat.shape[1])],
                                  *pL_vals)


    fig, ax = plt.subplots(4, 1, figsize=(8, 10), constrained_layout=True,
                              sharex=True)

    max_vel_S_z = np.max(np.abs(S_vel_N_np[:, 2]))
    max_pos_S_z = np.max(np.abs(S_pos_np[:, 2]))

    for i in (0, 1, 3, 4, 5, 6, 7, 8):
            ax[0].plot(sys.times, np.rad2deg(resultat[:, i]),
                       label=f'{bezeichnung[i]}')
            ax[1].plot(sys.times, resultat[:, i+9],
                       label=f'{bezeichnung[i+9]}')
            ax[1].axhline(0, color='k', lw=0.5, ls='--')

    ax[0].set_ylabel('[deg]')
    ax[0].set_title('Generalized Coordinates except $q_3, u_3$')
    ax[1].set_ylabel('[rad/s]')
    ax[1].set_title('Gen. speeds except $u_3$. Max. speed of S in N.z direction '
                f'is {max_vel_S_z:.2e} m/s')
    ax[-1].set_xlabel('Time [s]')
    ax[0].legend()
    ax[1].legend()

    koords = ['x', 'y', 'z']
    for i in range(3):
        ax[2].plot(sys.times, S_pos_np[:, i], label=f'S_pos_N_{koords[i]}')
    ax[2].set_ylabel('[m]')
    ax[2].set_title('Position of S in N. Max deviation of S from zero in '
                f'N.z direction is {max_pos_S_z:.2e} m')
    ax[2].set_xlabel('Time [s]')
    ax[2].legend()

    ax[3].plot(sys.times, np.rad2deg(resultat[:, 2]), label='$q_3$')
    ax[3].plot(sys.times, np.rad2deg(resultat[:, 11]), label='$u_3$')
    ax[3].set_title('Generalized coordinate $q_3, u_3$')
    ax[3].set_ylabel('deg / deg/sec')
    ax[3].legend()

    fig, ax = plt.subplots(figsize=(8, 8))
    min_x = np.min(S_pos_np[:, 0])
    max_x = np.max(S_pos_np[:, 0])
    min_y = np.min(S_pos_np[:, 1])
    max_y = np.max(S_pos_np[:, 1])
    ax.set_xlim(min_x - 0.01, max_x + 0.01)
    ax.set_ylim(min_y - 0.01, max_y + 0.01)
    ax.set_aspect('equal')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_aspect('equal')
    ax.set_title('Trajectory of S in the X/Y plane')

    points = np.array([S_pos_np[:, 0], S_pos_np[:, 1]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='inferno',
                        norm=plt.Normalize(sys.times[0], sys.times[-1]))
    lc.set_array(sys.times[:-1])
    lc.set_linewidth(0.75)
    ax.add_collection(lc)
    cbar = plt.colorbar(lc, ax=ax, shrink=0.7)
    cbar.set_label('Time [s]')


Plot Total Energy
-----------------

It drops as it should, since there is friction in the system.


.. jupyter-execute::

    pot_energy = m * g * (Ro.pos_from(O).dot(N.z))
    kin_energy = rattleback.kinetic_energy(N)
    pot_lam = sm.lambdify(qL1 + pL_pos, pot_energy, cse=True)
    kin_lam = sm.lambdify(qL1 + pL_pos, kin_energy, cse=True)

    pot_np = np.empty(resultat.shape[0])
    kin_np = np.empty(resultat.shape[0])
    total_np = np.empty(resultat.shape[0])

    for i in range(resultat.shape[0]):
        pot_np[i] = pot_lam(*resultat[i], *pL_vals)
        kin_np[i] = kin_lam(*resultat[i], *pL_vals)
        total_np[i] = pot_np[i] + kin_np[i]

    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 2), constrained_layout=True)
    ax2.plot(sys.times, total_np)
    ax2.set_ylabel('[Nm]')
    ax2.set_title('Total Energy. It drops, as it should, '
                  'due to friction in the system.')
    _ = ax2.set_xlabel('Time [s]')


Visualization
-------------

``groesse`` is an empirical value to control the size of the rattleback
in the animation.

.. jupyter-execute::

    groesse = 2.5

Define some points on the rattleback.

.. jupyter-execute::

    point1, point2, point3 = sm.symbols('point1, point2, point3', cls=me.Point)
    point1.set_pos(Oe, a/2 * R.x + c/2 * R.z)
    point2.set_pos(Oe, a * R.x)
    farben = ['red', 'green', 'blue']

    viz_frames = []

As the axes of pythreejs are fixed, and the equations of motion were set up
differently, they are rotated to match pythreejs' axes.

.. jupyter-execute::

    Bh = me.ReferenceFrame('Bh')
    Bh.orient_body_fixed(N, (sm.pi/2, 0, 0), 'XYZ')

Start the animation.

.. jupyter-execute::

    rattle_shape = Box(name='rattle',
                       width=2 * a,
                       height=2 * b,
                       depth=2 * c,
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

    scene = Scene(Bh, Oe, *viz_frames)

    scene.times = times
    sys.constants1 = {key: sys.constants[key] * groesse
                      for key in sys.constants.keys()}
    scene.constants = sys.constants1
    scene.states_symbols = q_ind + u_ind + u_dep
    scene.states_trajectories = resultat
    scene.display_jupyter(axes_arrow_length=5)
