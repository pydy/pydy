#

.. raw:: html

   <center>

Astrobee: A Holonomic Free-Flying Space Robot

.. raw:: html

   <center>

.. raw:: html

   <center>

.. raw:: html

   <video controls src="Astrobee.mov" width="500" />

.. raw:: html

   </center>

Dependencies
------------

.. code:: ipython3

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

.. code:: ipython3

    ISS = me.ReferenceFrame('N') # ISS RF
    B = me.ReferenceFrame('B') # body RF
    
    q1, q2, q3 = me.dynamicsymbols('q1:4') # attitude coordinates (Euler angles)
    
    B.orient(ISS, 'Body', (q1, q2, q3), 'xyz') # body RF

.. code:: ipython3

    t = me.dynamicsymbols._t

Significant Points
------------------

.. code:: ipython3

    O = me.Point('O') # fixed point in the ISS
    O.set_vel(ISS, 0)

.. code:: ipython3

    x, y, z = me.dynamicsymbols('x, y, z') # translation coordinates (position of the mass-center of Astrobee relative to 'O')
    l = sm.symbols('l') # length of Astrobee (side of cube)

.. code:: ipython3

    C = O.locatenew('C', x * ISS.x + y * ISS.y + z * ISS.z) # Astrobee CM

Kinematical Differential Equations
----------------------------------

.. code:: ipython3

    ux = me.dynamicsymbols('u_x')
    uy = me.dynamicsymbols('u_y')
    uz = me.dynamicsymbols('u_z')
    u1 = me.dynamicsymbols('u_1')
    u2 = me.dynamicsymbols('u_2')
    u3 = me.dynamicsymbols('u_3')

.. code:: ipython3

    z1 = sm.Eq(ux, x.diff())
    z2 = sm.Eq(uy, y.diff())
    z3 = sm.Eq(uz, z.diff())
    z4 = sm.Eq(u1, q1.diff())
    z5 = sm.Eq(u2, q2.diff())
    z6 = sm.Eq(u3, q3.diff())
    u = sm.solve([z1, z2, z3, z4, z5, z6], x.diff(), y.diff(), z.diff(), q1.diff(), q2.diff(), q3.diff())
    u




.. math::

    \displaystyle \left\{ \dot{q}_{1} : u_{1}, \  \dot{q}_{2} : u_{2}, \  \dot{q}_{3} : u_{3}, \  \dot{x} : u_{x}, \  \dot{y} : u_{y}, \  \dot{z} : u_{z}\right\}



Translational Motion
--------------------

Velocity
~~~~~~~~

.. code:: ipython3

    C.set_vel(ISS, C.pos_from(O).dt(ISS).subs(u))
    V_B_ISS_ISS = C.vel(ISS)
    V_B_ISS_ISS # "velocity of Astrobee CM w.r.t ISS RF expressed in ISS RF" 




.. math::

    \displaystyle u_{x}\mathbf{\hat{n}_x} + u_{y}\mathbf{\hat{n}_y} + u_{z}\mathbf{\hat{n}_z}



Acceleration
~~~~~~~~~~~~

.. code:: ipython3

    A_B_ISS_ISS = C.acc(ISS).subs(u) #.subs(ud)
    A_B_ISS_ISS # "acceleration of Astrobee CM w.r.t ISS RF expressed in ISS RF" 




.. math::

    \displaystyle \dot{u}_{x}\mathbf{\hat{n}_x} + \dot{u}_{y}\mathbf{\hat{n}_y} + \dot{u}_{z}\mathbf{\hat{n}_z}



Angular Motion
--------------

Angular Velocity
~~~~~~~~~~~~~~~~

.. code:: ipython3

    B.set_ang_vel(ISS, B.ang_vel_in(ISS).subs(u))
    Omega_B_ISS_B = B.ang_vel_in(ISS)
    Omega_B_ISS_B # "angular velocity of body RF w.r.t ISS RF expressed in body RF" 




.. math::

    \displaystyle (u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right))\mathbf{\hat{b}_x} + (- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right))\mathbf{\hat{b}_y} + (u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3})\mathbf{\hat{b}_z}



Angular Acceleration
~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    Alpha_B_ISS_B = B.ang_acc_in(ISS).subs(u) #.subs(ud)
    Alpha_B_ISS_B # "angular acceleration of body RF w.r.t ISS RF expressed in body RF" 




.. math::

    \displaystyle (- u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{1} u_{3} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} u_{3} \operatorname{cos}\left(q_{3}\right) + \operatorname{sin}\left(q_{3}\right) \dot{u}_{2} + \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) \dot{u}_{1})\mathbf{\hat{b}_x} + (u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{sin}\left(q_{3}\right) - u_{1} u_{3} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{2} u_{3} \operatorname{sin}\left(q_{3}\right) - \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) \dot{u}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{u}_{2})\mathbf{\hat{b}_y} + (u_{1} u_{2} \operatorname{cos}\left(q_{2}\right) + \operatorname{sin}\left(q_{2}\right) \dot{u}_{1} + \dot{u}_{3})\mathbf{\hat{b}_z}



Mass and Inertia
----------------

.. code:: ipython3

    m = sm.symbols('m') # Astrobee mass
    
    Ix, Iy, Iz = sm.symbols('I_x, I_y, I_z') # principal moments of inertia
    
    I = me.inertia(B, Ix, Iy, Iz) # inertia dyadic
    I




.. math::

    \displaystyle I_{x}\mathbf{\hat{b}_x}\otimes \mathbf{\hat{b}_x} + I_{y}\mathbf{\hat{b}_y}\otimes \mathbf{\hat{b}_y} + I_{z}\mathbf{\hat{b}_z}\otimes \mathbf{\hat{b}_z}



Loads
-----

Forces
~~~~~~

.. code:: ipython3

    Fx_mag, Fy_mag, Fz_mag = me.dynamicsymbols('Fmag_x, Fmag_y, Fmag_z')
    
    Fx = Fx_mag * ISS.x
    Fy = Fy_mag * ISS.y
    Fz = Fz_mag * ISS.z
    
    Fx, Fy, Fz




.. math::

    \displaystyle \left( \left|{F}\right|_{x}\mathbf{\hat{n}_x}, \  \left|{F}\right|_{y}\mathbf{\hat{n}_y}, \  \left|{F}\right|_{z}\mathbf{\hat{n}_z}\right)





.. code:: ipython3

    T1_mag, T2_mag, T3_mag = me.dynamicsymbols('Tmag_1, Tmag_2, Tmag_3')
    
    T1 = T1_mag * B.x
    T2 = T2_mag * B.y
    T3 = T3_mag * B.z
    
    T1, T2, T3




.. math::

    \displaystyle \left( \left|{T}\right|_{1}\mathbf{\hat{b}_x}, \  \left|{T}\right|_{2}\mathbf{\hat{b}_y}, \  \left|{T}\right|_{3}\mathbf{\hat{b}_z}\right)



Kane’s Method
-------------

.. code:: ipython3

    kdes = [z1.rhs - z1.lhs, z2.rhs - z2.lhs, z3.rhs - z3.lhs, z4.rhs - z4.lhs, z5.rhs - z5.lhs, z6.rhs - z6.lhs]

.. code:: ipython3

    body = me.RigidBody('body', C, B, m, (I, C))
    bodies = [body]

.. code:: ipython3

    loads = [
             (C, Fx),
             (C, Fy),
             (C, Fz),
             (B, T1),
             (B, T2),
             (B, T3)
            ]

.. code:: ipython3

    kane = me.KanesMethod(ISS, (x, y, z, q1, q2, q3), (ux, uy, uz, u1, u2, u3), kd_eqs=kdes)

.. code:: ipython3

    fr, frstar = kane.kanes_equations(bodies, loads=loads)

.. code:: ipython3

    # fr

.. code:: ipython3

    # frstar

Simulation
----------

.. code:: ipython3

    sys = System(kane)

.. code:: ipython3

    sys.constants_symbols




.. math::

    \displaystyle \left\{I_{x}, I_{y}, I_{z}, m\right\}



.. code:: ipython3

    sys.constants = {
                     Ix: 0.1083,
                     Iy: 0.1083,
                     Iz: 0.1083,
                     m: 7
                    }

.. code:: ipython3

    sys.constants




.. math::

    \displaystyle \left\{ I_{x} : 0.1083, \  I_{y} : 0.1083, \  I_{z} : 0.1083, \  m : 7\right\}



.. code:: ipython3

    sys.times = np.linspace(0.0, 50.0, num=1000)

.. code:: ipython3

    sys.coordinates




.. math::

    \displaystyle \left[ x, \  y, \  z, \  q_{1}, \  q_{2}, \  q_{3}\right]



.. code:: ipython3

    sys.speeds




.. math::

    \displaystyle \left[ u_{x}, \  u_{y}, \  u_{z}, \  u_{1}, \  u_{2}, \  u_{3}\right]



.. code:: ipython3

    sys.states




.. math::

    \displaystyle \left[ x, \  y, \  z, \  q_{1}, \  q_{2}, \  q_{3}, \  u_{x}, \  u_{y}, \  u_{z}, \  u_{1}, \  u_{2}, \  u_{3}\right]



.. code:: ipython3

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

.. code:: ipython3

    sys.specifieds_symbols




.. math::

    \displaystyle \left\{\left|{F}\right|_{x}, \left|{F}\right|_{y}, \left|{F}\right|_{z}, \left|{T}\right|_{1}, \left|{T}\right|_{2}, \left|{T}\right|_{3}\right\}



.. code:: ipython3

    sys.specifieds = {
                      Fx_mag: 0.0,
                      Fy_mag: 0.0,
                      Fz_mag: 0.0,
                      T1_mag: 0.0,
                      T2_mag: 0.0,
                      T3_mag: 0.0
                     }

.. code:: ipython3

    states = sys.integrate()

.. code:: ipython3

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

.. code:: ipython3

    fig, ax = plt.subplots()
    ax.plot(sys.times, states)
    ax.set_xlabel('{} [s]'.format(sm.latex(t, mode='inline')));
    ax.set_ylabel('States');
    ax.legend(['$x$', '$y$', '$z$', '$q_1$', '$q_2$', '$q_3$', '$u_x$', '$u_y$', '$u_z$', '$u_1$', '$u_2$', '$u_3$'], fontsize=10)
    plt.show()



.. image:: Space_Robot_files/Space_Robot_53_0.png


3D Visualization
----------------

.. code:: ipython3

    from pydy.viz.shapes import Cube, Cylinder, Sphere, Plane
    from pydy.viz.visualization_frame import VisualizationFrame
    from pydy.viz import Scene
    from ipywidgets import Image, Video
    import pythreejs as pjs
    from stl import mesh

.. code:: ipython3

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

Linearization
-------------

.. code:: ipython3

    f = fr + frstar
    f




.. math::

    \displaystyle \left[\begin{matrix}- m \dot{u}_{x} + \left|{F}\right|_{x}\\- m \dot{u}_{y} + \left|{F}\right|_{y}\\- m \dot{u}_{z} + \left|{F}\right|_{z}\\- I_{z} \operatorname{sin}\left(q_{2}\right) \dot{u}_{3} - \left(I_{x} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - I_{y} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{u}_{2} - \left(I_{x} \left(- u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{1} u_{3} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} u_{3} \operatorname{cos}\left(q_{3}\right)\right) - I_{y} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) + I_{z} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right)\right) \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + \left(I_{x} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) + I_{y} \left(u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{sin}\left(q_{3}\right) - u_{1} u_{3} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{2} u_{3} \operatorname{sin}\left(q_{3}\right)\right) - I_{z} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right)\right) \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) - \left(- I_{x} \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) + I_{y} \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) + I_{z} u_{1} u_{2} \operatorname{cos}\left(q_{2}\right)\right) \operatorname{sin}\left(q_{2}\right) - \left(I_{x} \operatorname{cos}^{2}\left(q_{2}\right) \operatorname{cos}^{2}\left(q_{3}\right) + I_{y} \operatorname{sin}^{2}\left(q_{3}\right) \operatorname{cos}^{2}\left(q_{2}\right) + I_{z} \operatorname{sin}^{2}\left(q_{2}\right)\right) \dot{u}_{1} + \left|{T}\right|_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - \left|{T}\right|_{2} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + \left|{T}\right|_{3} \operatorname{sin}\left(q_{2}\right)\\- \left(I_{x} \operatorname{sin}^{2}\left(q_{3}\right) + I_{y} \operatorname{cos}^{2}\left(q_{3}\right)\right) \dot{u}_{2} - \left(I_{x} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - I_{y} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{u}_{1} - \left(I_{x} \left(- u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{1} u_{3} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} u_{3} \operatorname{cos}\left(q_{3}\right)\right) - I_{y} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) + I_{z} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right)\right) \operatorname{sin}\left(q_{3}\right) - \left(I_{x} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) + I_{y} \left(u_{1} u_{2} \operatorname{sin}\left(q_{2}\right) \operatorname{sin}\left(q_{3}\right) - u_{1} u_{3} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - u_{2} u_{3} \operatorname{sin}\left(q_{3}\right)\right) - I_{z} \left(u_{1} \operatorname{sin}\left(q_{2}\right) + u_{3}\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right)\right) \operatorname{cos}\left(q_{3}\right) + \left|{T}\right|_{1} \operatorname{sin}\left(q_{3}\right) + \left|{T}\right|_{2} \operatorname{cos}\left(q_{3}\right)\\I_{x} \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) - I_{y} \left(- u_{1} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + u_{2} \operatorname{cos}\left(q_{3}\right)\right) \left(u_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) + u_{2} \operatorname{sin}\left(q_{3}\right)\right) - I_{z} u_{1} u_{2} \operatorname{cos}\left(q_{2}\right) - I_{z} \operatorname{sin}\left(q_{2}\right) \dot{u}_{1} - I_{z} \dot{u}_{3} + \left|{T}\right|_{3}\end{matrix}\right]



.. code:: ipython3

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

.. code:: ipython3

    f_lin = f.subs(V) + f.jacobian(V_keys).subs(V)*(V_keys - V_values)

.. code:: ipython3

    # sm.simplify(f)

.. code:: ipython3

    sm.simplify(f.subs(sys.constants))




.. math::

    \displaystyle \left[\begin{matrix}\left|{F}\right|_{x} - 7 \dot{u}_{x}\\\left|{F}\right|_{y} - 7 \dot{u}_{y}\\\left|{F}\right|_{z} - 7 \dot{u}_{z}\\1.0 \left|{T}\right|_{1} \operatorname{cos}\left(q_{2}\right) \operatorname{cos}\left(q_{3}\right) - 1.0 \left|{T}\right|_{2} \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{2}\right) + 1.0 \left|{T}\right|_{3} \operatorname{sin}\left(q_{2}\right) - 0.1083 u_{2} u_{3} \operatorname{cos}\left(q_{2}\right) - 0.1083 \operatorname{sin}\left(q_{2}\right) \dot{u}_{3} - 0.1083 \dot{u}_{1}\\1.0 \left|{T}\right|_{1} \operatorname{sin}\left(q_{3}\right) + 1.0 \left|{T}\right|_{2} \operatorname{cos}\left(q_{3}\right) + 0.1083 u_{1} u_{3} \operatorname{cos}\left(q_{2}\right) - 0.1083 \dot{u}_{2}\\\left|{T}\right|_{3} - 0.1083 u_{1} u_{2} \operatorname{cos}\left(q_{2}\right) - 0.1083 \operatorname{sin}\left(q_{2}\right) \dot{u}_{1} - 0.1083 \dot{u}_{3}\end{matrix}\right]



.. code:: ipython3

    us = sm.Matrix([ux, uy, uz, u1, u2, u3])
    us_diff = sm.Matrix([ux.diff(), uy.diff(), uz.diff(), u1.diff(), u2.diff(), u3.diff()])
    qs = sm.Matrix([x, y, z, q1, q2, q3])
    rs = sm.Matrix([Fx_mag, Fy_mag, Fz_mag, T1_mag, T2_mag, T3_mag])

If :math:`f_{lin}` is used, :math:`M_l \rightarrow` singular

:math:`\because` inversion of :math:`M_l` is required, use :math:`f` and
then substitute for :math:`V`

.. code:: ipython3

    Ml = f.jacobian(us_diff).subs(sys.constants).subs(V)
    Ml




.. math::

    \displaystyle \left[\begin{matrix}-7 & 0 & 0 & 0 & 0 & 0\\0 & -7 & 0 & 0 & 0 & 0\\0 & 0 & -7 & 0 & 0 & 0\\0 & 0 & 0 & -0.1083 & 0 & 0\\0 & 0 & 0 & 0 & -0.1083 & 0\\0 & 0 & 0 & 0 & 0 & -0.1083\end{matrix}\right]



.. code:: ipython3

    Cl = f.jacobian(us).subs(V)
    Cl.subs(sys.constants)




.. math::

    \displaystyle \left[\begin{matrix}0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]



.. code:: ipython3

    Kl = f.jacobian(qs).subs(V)
    sm.simplify(Kl.subs(sys.constants))




.. math::

    \displaystyle \left[\begin{matrix}0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]



.. code:: ipython3

    Hl = -f.jacobian(rs).subs(V)
    sm.simplify(Hl.subs(sys.constants))




.. math::

    \displaystyle \left[\begin{matrix}-1 & 0 & 0 & 0 & 0 & 0\\0 & -1 & 0 & 0 & 0 & 0\\0 & 0 & -1 & 0 & 0 & 0\\0 & 0 & 0 & -1 & 0 & 0\\0 & 0 & 0 & 0 & -1 & 0\\0 & 0 & 0 & 0 & 0 & -1\end{matrix}\right]



.. code:: ipython3

    A = sm.Matrix([[(-Ml.inv()*Cl), (-Ml.inv()*Kl)], [(sm.eye(6)), sm.zeros(6, 6)]])
    sm.simplify(A.subs(sys.constants))




.. math::

    \displaystyle \left[\begin{array}{cccccccccccc}0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\end{array}\right]



.. code:: ipython3

    sm.simplify(A).subs(sys.constants)*(us.col_join(qs))




.. math::

    \displaystyle \left[\begin{matrix}0\\0\\0\\0\\0\\0\\u_{x}\\u_{y}\\u_{z}\\u_{1}\\u_{2}\\u_{3}\end{matrix}\right]



.. code:: ipython3

    B = sm.Matrix([[Ml.inv() * Hl], [sm.zeros(6, 6)]])
    sm.nsimplify(B.subs(sys.constants))




.. math::

    \displaystyle \left[\begin{matrix}\frac{1}{7} & 0 & 0 & 0 & 0 & 0\\0 & \frac{1}{7} & 0 & 0 & 0 & 0\\0 & 0 & \frac{1}{7} & 0 & 0 & 0\\0 & 0 & 0 & 9.23361034164358 & 0 & 0\\0 & 0 & 0 & 0 & 9.23361034164358 & 0\\0 & 0 & 0 & 0 & 0 & 9.23361034164358\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]



.. code:: ipython3

    sm.simplify(B).subs(sys.constants)*(rs)




.. math::

    \displaystyle \left[\begin{matrix}\frac{\left|{F}\right|_{x}}{7}\\\frac{\left|{F}\right|_{y}}{7}\\\frac{\left|{F}\right|_{z}}{7}\\9.23361034164358 \left|{T}\right|_{1}\\9.23361034164358 \left|{T}\right|_{2}\\9.23361034164358 \left|{T}\right|_{3}\\0\\0\\0\\0\\0\\0\end{matrix}\right]



.. code:: ipython3

    us.col_join(qs), (sm.simplify(A).subs(sys.constants)*(us.col_join(qs)) + sm.simplify(B).subs(sys.constants)*(rs)) # (x, Ax + Bu) => x_dot = Ax + Bu?




.. math::

    \displaystyle \left( \left[\begin{matrix}u_{x}\\u_{y}\\u_{z}\\u_{1}\\u_{2}\\u_{3}\\x\\y\\z\\q_{1}\\q_{2}\\q_{3}\end{matrix}\right], \  \left[\begin{matrix}\frac{\left|{F}\right|_{x}}{7}\\\frac{\left|{F}\right|_{y}}{7}\\\frac{\left|{F}\right|_{z}}{7}\\9.23361034164358 \left|{T}\right|_{1}\\9.23361034164358 \left|{T}\right|_{2}\\9.23361034164358 \left|{T}\right|_{3}\\u_{x}\\u_{y}\\u_{z}\\u_{1}\\u_{2}\\u_{3}\end{matrix}\right]\right)



.. code:: ipython3

    (us.col_join(qs))




.. math::

    \displaystyle \left[\begin{matrix}u_{x}\\u_{y}\\u_{z}\\u_{1}\\u_{2}\\u_{3}\\x\\y\\z\\q_{1}\\q_{2}\\q_{3}\end{matrix}\right]



.. code:: ipython3

    renderer



.. parsed-literal::

    Renderer(camera=PerspectiveCamera(aspect=0.75, position=(1.0, 1.0, 1.0), quaternion=(0.0, 0.0, 0.0, 1.0), scal…


.. code:: ipython3

    clip = pjs.AnimationClip(tracks=[body_m_track, body_l_track, body_r_track], duration=sys.times[-1])
    
    
    action = pjs.AnimationAction(pjs.AnimationMixer(scene_pjs), clip, scene_pjs)
    action



.. parsed-literal::

    AnimationAction(clip=AnimationClip(duration=50.0, tracks=(VectorKeyframeTrack(name='scene/Body_m.matrix', time…

