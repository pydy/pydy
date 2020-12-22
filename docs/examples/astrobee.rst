=============================================
Astrobee: A Holonomic Free-Flying Space Robot
=============================================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`astrobee` or Jupyter notebook:
   :jupyter-download:notebook:`astrobee`.

Astrobee is a new generation of free-flying robots aboard the International
Space Station (ISS). It is a cubic robot with sides measuring about 30 cm each.
The robot is propelled by two fans located on the sides of the robot and
servo-actuated louvred vent nozzles, which allow for full six degree-of-freedom
holonomic control [Smith2016]_. Here, the nonlinear dynamics of Astrobee are
modeled using Kane's method and the holonomic behavior of the system is
demonstrated. After derivation of the nonlinear equations of motion, the system
is linearized about a chosen operating point to obtain an explicit first order
state-space representation, which can be used for control design.

.. jupyter-execute::

    import sympy as sm
    import sympy.physics.mechanics as me
    from pydy.system import System
    import numpy as np
    import matplotlib.pyplot as plt
    from pydy.codegen.ode_function_generators import generate_ode_function
    from scipy.integrate import odeint
    import scipy.io as sio
    me.init_vprinting()

Reference Frames
----------------

.. jupyter-execute::

    ISS = me.ReferenceFrame('N') # ISS RF
    B = me.ReferenceFrame('B') # body RF

    q1, q2, q3 = me.dynamicsymbols('q1:4') # attitude coordinates (Euler angles)

    B.orient(ISS, 'Body', (q1, q2, q3), 'xyz') # body RF

    t = me.dynamicsymbols._t

Significant Points
------------------

.. jupyter-execute::

    O = me.Point('O') # fixed point in the ISS
    O.set_vel(ISS, 0)

    x, y, z = me.dynamicsymbols('x, y, z') # translation coordinates (position of the mass-center of Astrobee relative to 'O')
    l = sm.symbols('l') # length of Astrobee (side of cube)

    C = O.locatenew('C', x * ISS.x + y * ISS.y + z * ISS.z) # Astrobee CM

Kinematical Differential Equations
----------------------------------

.. jupyter-execute::

    ux = me.dynamicsymbols('u_x')
    uy = me.dynamicsymbols('u_y')
    uz = me.dynamicsymbols('u_z')
    u1 = me.dynamicsymbols('u_1')
    u2 = me.dynamicsymbols('u_2')
    u3 = me.dynamicsymbols('u_3')

    z1 = sm.Eq(ux, x.diff())
    z2 = sm.Eq(uy, y.diff())
    z3 = sm.Eq(uz, z.diff())
    z4 = sm.Eq(u1, q1.diff())
    z5 = sm.Eq(u2, q2.diff())
    z6 = sm.Eq(u3, q3.diff())
    u = sm.solve([z1, z2, z3, z4, z5, z6], x.diff(), y.diff(), z.diff(), q1.diff(), q2.diff(), q3.diff())
    u

Translational Motion
--------------------

Velocity
~~~~~~~~

.. jupyter-execute::

    C.set_vel(ISS, C.pos_from(O).dt(ISS).subs(u))
    V_B_ISS_ISS = C.vel(ISS)
    V_B_ISS_ISS # "velocity of Astrobee CM w.r.t ISS RF expressed in ISS RF"

Acceleration
~~~~~~~~~~~~

.. jupyter-execute::

    A_B_ISS_ISS = C.acc(ISS).subs(u) #.subs(ud)
    A_B_ISS_ISS # "acceleration of Astrobee CM w.r.t ISS RF expressed in ISS RF"

Angular Motion
--------------

Angular Velocity
~~~~~~~~~~~~~~~~

.. jupyter-execute::

    B.set_ang_vel(ISS, B.ang_vel_in(ISS).subs(u))
    Omega_B_ISS_B = B.ang_vel_in(ISS)
    Omega_B_ISS_B # "angular velocity of body RF w.r.t ISS RF expressed in body RF"

Angular Acceleration
~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    Alpha_B_ISS_B = B.ang_acc_in(ISS).subs(u) #.subs(ud)
    Alpha_B_ISS_B # "angular acceleration of body RF w.r.t ISS RF expressed in body RF"

Mass and Inertia
----------------

.. jupyter-execute::

    m = sm.symbols('m') # Astrobee mass

    Ix, Iy, Iz = sm.symbols('I_x, I_y, I_z') # principal moments of inertia

    I = me.inertia(B, Ix, Iy, Iz) # inertia dyadic
    I

Loads
-----

Forces
~~~~~~

.. jupyter-execute::

    Fx_mag, Fy_mag, Fz_mag = me.dynamicsymbols('Fmag_x, Fmag_y, Fmag_z')

    Fx = Fx_mag * ISS.x
    Fy = Fy_mag * ISS.y
    Fz = Fz_mag * ISS.z

    Fx, Fy, Fz

Torques
~~~~~~~

.. jupyter-execute::

    T1_mag, T2_mag, T3_mag = me.dynamicsymbols('Tmag_1, Tmag_2, Tmag_3')

    T1 = T1_mag * B.x
    T2 = T2_mag * B.y
    T3 = T3_mag * B.z

    T1, T2, T3

Kane’s Method
-------------

.. jupyter-execute::

    kdes = [z1.rhs - z1.lhs,
            z2.rhs - z2.lhs,
            z3.rhs - z3.lhs,
            z4.rhs - z4.lhs,
            z5.rhs - z5.lhs,
            z6.rhs - z6.lhs]

    body = me.RigidBody('body', C, B, m, (I, C))
    bodies = [body]

    loads = [
             (C, Fx),
             (C, Fy),
             (C, Fz),
             (B, T1),
             (B, T2),
             (B, T3)
            ]

    kane = me.KanesMethod(ISS, (x, y, z, q1, q2, q3), (ux, uy, uz, u1, u2, u3), kd_eqs=kdes)

    fr, frstar = kane.kanes_equations(bodies, loads=loads)

Simulation
----------

.. jupyter-execute::

    sys = System(kane)

    sys.constants_symbols

    sys.constants = {
                     Ix: 0.1083,
                     Iy: 0.1083,
                     Iz: 0.1083,
                     m: 7
                    }

    sys.constants

.. jupyter-execute::

    sys.times = np.linspace(0.0, 50.0, num=1000)

    sys.coordinates

.. jupyter-execute::

    sys.speeds

.. jupyter-execute::

    sys.states

.. jupyter-execute::

    sys.initial_conditions = {
                              x: 0.0,
                              y: 0.0,
                              z: 0.0,
                              q1: 0.0,
                              q2: 0.0,
                              q3: 0.0,
                              ux: 0.2,
                              uy: 0.0,
                              uz: 0.0,
                              u1: 0.0,
                              u2: 0.0,
                              u3: 0.5
                             }

.. jupyter-execute::

    sys.specifieds_symbols


.. jupyter-execute::

    sys.specifieds = {
                      Fx_mag: 0.0,
                      Fy_mag: 0.0,
                      Fz_mag: 0.0,
                      T1_mag: 0.0,
                      T2_mag: 0.0,
                      T3_mag: 0.0
                     }

.. jupyter-execute::

    states = sys.integrate()

.. jupyter-execute::

    import matplotlib.pyplot as plt

.. jupyter-execute::

    fig, ax = plt.subplots()
    ax.plot(sys.times, states)
    ax.set_xlabel('{} [s]'.format(sm.latex(t, mode='inline')));
    ax.set_ylabel('States');
    ax.legend(['$x$', '$y$', '$z$', '$q_1$', '$q_2$', '$q_3$', '$u_x$', '$u_y$', '$u_z$', '$u_1$', '$u_2$', '$u_3$'], fontsize=10)
    plt.show()

3D Visualization
----------------

.. jupyter-execute::

    from pydy.viz import Box, Cube, Sphere, Cylinder, VisualizationFrame, Scene

.. jupyter-execute::

    l = 0.32

    body_m_shape = Box(l, (1/2) * l, (2/3) * l, color='black', name='body_m_shape')
    body_l_shape = Box(l, (1/4) * l, l, color='green', name='body_l_shape')
    body_r_shape = Box(l, (1/4) * l, l, color='green', name='body_r_shape')

    v1 = VisualizationFrame('Body_m',
                            B,
                            C.locatenew('C_m', (1/6) * l * B.z),
                            body_m_shape)

    v2 = VisualizationFrame('Body_l',
                            B,
                            C.locatenew('C_l', (3/8) * l * -B.y),
                            body_l_shape)

    v3 = VisualizationFrame('Body_r',
                            B,
                            C.locatenew('C_r', (3/8) * l * B.y),
                            body_r_shape)

    scene = Scene(ISS, O, system=sys)

    scene.visualization_frames = [v1, v2, v3]


.. jupyter-execute::

   scene.display_jupyter(axes_arrow_length=1.0)

Linearization
-------------

.. jupyter-execute::

    f = fr + frstar
    f

.. jupyter-execute::

    V = {
          x: 0.0,
          y: 0.0,
          z: 0.0,
          q1: 0.0,
          q2: 0.0,
          q3: 0.0,
          ux: 0.0,
          uy: 0.0,
          uz: 0.0,
          u1: 0.0,
          u2: 0.0,
          u3: 0.0,
          Fx_mag: 0.0,
          Fy_mag: 0.0,
          Fz_mag: 0.0,
          T1_mag: 0.0,
          T2_mag: 0.0,
          T3_mag: 0.0
    }

    V_keys = sm.Matrix([ v for v in V.keys() ])
    V_values = sm.Matrix([ v for v in V.values() ])

.. jupyter-execute::

    us = sm.Matrix([ux, uy, uz, u1, u2, u3])
    us_diff = sm.Matrix([ux.diff(), uy.diff(), uz.diff(), u1.diff(), u2.diff(), u3.diff()])
    qs = sm.Matrix([x, y, z, q1, q2, q3])
    rs = sm.Matrix([Fx_mag, Fy_mag, Fz_mag, T1_mag, T2_mag, T3_mag])

.. jupyter-execute::

    Ml = f.jacobian(us_diff).subs(sys.constants).subs(V)
    Ml

.. jupyter-execute::

    Cl = f.jacobian(us).subs(V)
    Cl.subs(sys.constants)

.. jupyter-execute::

    Kl = f.jacobian(qs).subs(V)
    sm.simplify(Kl.subs(sys.constants))

.. jupyter-execute::

    Hl = -f.jacobian(rs).subs(V)
    sm.simplify(Hl.subs(sys.constants))

.. jupyter-execute::

    A = sm.Matrix([[(-Ml.inv()*Cl), (-Ml.inv()*Kl)], [(sm.eye(6)), sm.zeros(6, 6)]])
    sm.simplify(A.subs(sys.constants))

.. jupyter-execute::

    B = sm.Matrix([[Ml.inv() * Hl], [sm.zeros(6, 6)]])
    sm.nsimplify(B.subs(sys.constants))

References
----------

.. [Smith2016] Smith, T., Barlow, J., Bualat, M., Fong, T., Provencher, C.,
   Sanchez, H., & Smith, E. (2016). Astrobee: A new platform for free-flying
   robotics on the international space station.
