=============================================
Astrobee: A Holonomic Free-Flying Space Robot
=============================================

[Bualat2015]_

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

.. jupyter-execute::

    t = me.dynamicsymbols._t

Significant Points
------------------

.. jupyter-execute::

    O = me.Point('O') # fixed point in the ISS
    O.set_vel(ISS, 0)

.. jupyter-execute::

    x, y, z = me.dynamicsymbols('x, y, z') # translation coordinates (position of the mass-center of Astrobee relative to 'O')
    l = sm.symbols('l') # length of Astrobee (side of cube)

.. jupyter-execute::

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

.. jupyter-execute::

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

    kdes = [z1.rhs - z1.lhs, z2.rhs - z2.lhs, z3.rhs - z3.lhs, z4.rhs - z4.lhs, z5.rhs - z5.lhs, z6.rhs - z6.lhs]

.. jupyter-execute::

    body = me.RigidBody('body', C, B, m, (I, C))
    bodies = [body]

.. jupyter-execute::

    loads = [
             (C, Fx),
             (C, Fy),
             (C, Fz),
             (B, T1),
             (B, T2),
             (B, T3)
            ]

.. jupyter-execute::

    kane = me.KanesMethod(ISS, (x, y, z, q1, q2, q3), (ux, uy, uz, u1, u2, u3), kd_eqs=kdes)

.. jupyter-execute::

    fr, frstar = kane.kanes_equations(bodies, loads=loads)



Simulation
----------

.. jupyter-execute::

    sys = System(kane)

.. jupyter-execute::

    sys.constants_symbols



.. jupyter-execute::

    sys.constants = {
                     Ix: 0.1083,
                     Iy: 0.1083,
                     Iz: 0.1083,
                     m: 7
                    }

.. jupyter-execute::

    sys.constants



.. jupyter-execute::

    sys.times = np.linspace(0.0, 50.0, num=1000)

.. jupyter-execute::

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

    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 200
    mpl.rc('font',**{'family':'serif','sans-serif':['Computer Modern Roman']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    mpl.rc('text', usetex=True)
    from matplotlib.pyplot import cm
    color=cm.rainbow(np.linspace(0,1,12))
    from cycler import cycler
    mpl.rcParams['axes.prop_cycle'] = cycler(color=color)
    mpl.rcParams.update({'figure.autolayout': True})
    mpl.rcParams.update({'font.size': 12})
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D, art3d

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

    from pydy.viz.shapes import Cube, Cylinder, Sphere, Plane
    from pydy.viz.visualization_frame import VisualizationFrame
    from pydy.viz import Scene
    from ipywidgets import Image, Video
    import pythreejs as pjs
    from stl import mesh

.. jupyter-execute::

    l = 0.32
    
    body_m_shape = Cube(l, color='black')
    body_l_shape = Cube(l, color='green')
    body_r_shape = Cube(l, color='green')
    
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
    
    scene = Scene(ISS, O, v1, v2, v3, system=sys)
    scene.create_static_html(overwrite=True, silent=True)
    
    body_m_mesh = pjs.Mesh(
        pjs.BoxBufferGeometry(l, (1/2) * l, (2/3) * l),
        pjs.MeshStandardMaterial(color='black'),
        name="Body_m"
    )
    
    body_l_mesh = pjs.Mesh(
        pjs.BoxBufferGeometry(l, (1/4) * l, l),
        pjs.MeshStandardMaterial(color='green'),
        name="Body_l"
    )
    
    body_r_mesh = pjs.Mesh(
        pjs.BoxBufferGeometry(l, (1/4) * l, l),
        pjs.MeshStandardMaterial(color='green'),
        name="Body_r"
    )
    
    body_m_matrices = v1.evaluate_transformation_matrix(states, list(sys.constants.values()))
    body_l_matrices = v2.evaluate_transformation_matrix(states, list(sys.constants.values()))
    body_r_matrices = v3.evaluate_transformation_matrix(states, list(sys.constants.values()))
    
    body_m_track = pjs.VectorKeyframeTrack(
        name='scene/Body_m.matrix',
        times=list(sys.times),
        values=body_m_matrices)
    
    body_l_track = pjs.VectorKeyframeTrack(
        name='scene/Body_l.matrix',
        times=list(sys.times),
        values=body_l_matrices)
    
    body_r_track = pjs.VectorKeyframeTrack(
        name='scene/Body_r.matrix',
        times=list(sys.times),
        values=body_r_matrices)
    
    body_m_mesh.matrixAutoUpdate = False
    body_l_mesh.matrixAutoUpdate = False
    body_r_mesh.matrixAutoUpdate = False
    
    body_m_mesh.matrix = body_m_matrices[0]
    body_l_mesh.matrix = body_l_matrices[0]
    body_r_mesh.matrix = body_r_matrices[0]
    
    x_arrow = pjs.ArrowHelper(dir=[1, 0, 0], length=0.75, color='blue')
    y_arrow = pjs.ArrowHelper(dir=[0, 1, 0], length=0.75, color='red')
    z_arrow = pjs.ArrowHelper(dir=[0, 0, 1], length=0.75,color='green')
    
    view_width = 960
    view_height = 720
    
    camera = pjs.PerspectiveCamera(position=[1, 1, 1],
                                   aspect=view_width/view_height)
    key_light = pjs.DirectionalLight(position=[1, 1, 0])
    ambient_light = pjs.AmbientLight()
    
    scene_pjs = pjs.Scene(children=[body_m_mesh, body_l_mesh, body_r_mesh,
                                    x_arrow, y_arrow, z_arrow, 
                                    camera, key_light, ambient_light])
    
    controller = pjs.OrbitControls(controlling=camera)
    renderer = pjs.Renderer(camera=camera, scene=scene_pjs, controls=[controller], width=view_width, height=view_height)


.. jupyter-execute::

    renderer


.. jupyter-execute::

    clip = pjs.AnimationClip(tracks=[body_m_track, body_l_track, body_r_track], duration=sys.times[-1])
    
    
    action = pjs.AnimationAction(pjs.AnimationMixer(scene_pjs), clip, scene_pjs)
    action



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

.. [Bualat2015] Bualat, M., Barlow, J., Fong, T., Provencher, C. and Smith, T., 2015. Astrobee: Developing a free-flying robot for the international space station. In AIAA SPACE 2015 Conference and Exposition (p. 4643).
