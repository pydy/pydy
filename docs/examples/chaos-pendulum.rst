==============
Chaos Pendulum
==============

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`chaos-pendulum` or Jupyter notebook:
   :jupyter-download:notebook:`chaos-pendulum`.

This example gives a simple demostration of chaotic behavior in a simple
two body system. The system is made up of a slender rod that is
connected to the ceiling at one end with a revolute joint that rotates
about the :math:`\hat{\mathbf{n}}_y` unit vector. At the other end of
the rod a flat plate is attached via a second revolute joint allowing
the plate to rotate about the rod’s axis with aligns with the
:math:`\hat{\mathbf{a}_z}` unit vector.

.. image:: chaos-pendulum.svg

Setup
=====

.. jupyter-execute::

    import numpy as np
    import matplotlib.pyplot as plt
    import sympy as sm
    import sympy.physics.mechanics as me
    from pydy.system import System
    from pydy.viz import Cylinder, Plane, VisualizationFrame, Scene

.. jupyter-execute::

    %matplotlib inline

.. jupyter-execute::

    me.init_vprinting(use_latex='mathjax')

Define Variables
================

First define the system constants:

-  :math:`m_A`: Mass of the slender rod.
-  :math:`m_B`: Mass of the plate.
-  :math:`l_B`: Distance from :math:`N_o` to :math:`B_o` along the
   slender rod’s axis.
-  :math:`w`: The width of the plate.
-  :math:`h`: The height of the plate.
-  :math:`g`: The acceleratoin due to gravity.

.. jupyter-execute::

    mA, mB, lB, w, h, g = sm.symbols('m_A, m_B, L_B, w, h, g')

There are two time varying generalized coordinates:

-  :math:`\theta(t)`: The angle of the slender rod with respect to the
   ceiling.
-  :math:`\phi(t)`: The angle of the plate with respect to the slender
   rod.

The two generalized speeds will then be defined as:

-  :math:`\omega(t)=\dot{\theta}`: The angular rate of the slender rod
   with respect to the ceiling.
-  :math:`\alpha(t)=\dot{\phi}`: The angluer rate of the plate with
   respect to the slender rod.

.. jupyter-execute::

    theta, phi = me.dynamicsymbols('theta, phi')
    omega, alpha = me.dynamicsymbols('omega, alpha')

The kinematical differential equations are defined in this fashion for
the ``KanesMethod`` class:

.. math::

   0 = \omega - \dot{\theta}\\
   0 = \alpha - \dot{\phi}

.. jupyter-execute::

    kin_diff = (omega - theta.diff(), alpha - phi.diff())
    kin_diff

Define Orientations
===================

There are three reference frames. These are defined as such:

.. jupyter-execute::

    N = me.ReferenceFrame('N')
    A = me.ReferenceFrame('A')
    B = me.ReferenceFrame('B')

The frames are oriented with respect to each other by simple revolute
rotations. The following lines set the orientations:

.. jupyter-execute::

    A.orient(N, 'Axis', (theta, N.y))
    B.orient(A, 'Axis', (phi, A.z))

Define Positions
================

Three points are necessary to define the problem:

-  :math:`N_o`: The fixed point which the slender rod rotates about.
-  :math:`A_o`: The center of mass of the slender rod.
-  :math:`B_o`: The center of mass of the plate.

.. jupyter-execute::

    No = me.Point('No')
    Ao = me.Point('Ao')
    Bo = me.Point('Bo')

The two centers of mass positions can be set relative to the fixed
point, :math:`N_o`.

.. jupyter-execute::

    lA = (lB - h / 2) / 2
    Ao.set_pos(No, lA * A.z)
    Bo.set_pos(No, lB * A.z)

Specify the Velocities
======================

The generalized speeds should be used in the definition of the linear
and angular velocities when using Kane’s method. For simple rotations
and the defined kinematical differential equations the angular rates
are:

.. jupyter-execute::

    A.set_ang_vel(N, omega * N.y)
    B.set_ang_vel(A, alpha * A.z)

Once the angular velocities are specified the linear velocities can be
computed using the two point velocity thereom, starting with the origin
point having a velocity of zero.

.. jupyter-execute::

    No.set_vel(N, 0)

.. jupyter-execute::

    Ao.v2pt_theory(No, N, A)

.. jupyter-execute::

    Bo.v2pt_theory(No, N, A)

Inertia
=======

The central inertia of the symmetric slender rod with respect to its
reference frame is a function of its length and its mass.

.. jupyter-execute::

    IAxx = sm.S(1) / 12 * mA * (2 * lA)**2
    IAyy = IAxx
    IAzz = 0

    IA = (me.inertia(A, IAxx, IAyy, IAzz), Ao)

This gives the inertia tensor:

.. jupyter-execute::

    IA[0].to_matrix(A)

The central inerita of the symmetric plate with respect to its reference
frame is a function of its width and height.

.. jupyter-execute::

    IBxx = sm.S(1)/12 * mB * h**2
    IByy = sm.S(1)/12 * mB * (w**2 + h**2)
    IBzz = sm.S(1)/12 * mB * w**2

    IB = (me.inertia(B, IBxx, IByy, IBzz), Bo)

.. jupyter-execute::

    IB[0].to_matrix(B)

All of the information to define the two rigid bodies are now available.
This information is used to create an object for the rod and the plate.

.. jupyter-execute::

    rod = me.RigidBody('rod', Ao, A, mA, IA)

.. jupyter-execute::

    plate = me.RigidBody('plate', Bo, B, mB, IB)

Loads
=====

The only loads in this problem is the force due to gravity that acts on
the center of mass of each body. These forces are specified with a tuple
containing the point of application and the force vector.

.. jupyter-execute::

    rod_gravity = (Ao, mA * g * N.z)
    plate_gravity = (Bo, mB * g * N.z)

Equations of motion
===================

Now that the kinematics, kinetics, and inertia have all been defined the
``KanesMethod`` class can be used to generate the equations of motion of
the system. In this case the independent generalized speeds, independent
generalized speeds, the kinematical differential equations, and the
inertial reference frame are used to initialize the class.

.. jupyter-execute::

    kane = me.KanesMethod(N, q_ind=(theta, phi), u_ind=(omega, alpha), kd_eqs=kin_diff)

The equations of motion are then generated by passing in all of the
loads and bodies to the ``kanes_equations`` method. This produces
:math:`f_r` and :math:`f_r^*`.

.. jupyter-execute::

    bodies = (rod, plate)
    loads = (rod_gravity, plate_gravity)

    fr, frstar = kane.kanes_equations(bodies, loads)

.. jupyter-execute::

    sm.trigsimp(fr)

.. jupyter-execute::

    sm.trigsimp(frstar)

Simulation
==========

The equations of motion can now be simulated numerically. Values for the
constants, initial conditions, and time are provided to the ``System``
class along with the symbolic ``KanesMethod`` object.

.. jupyter-execute::

    sys = System(kane)

.. jupyter-execute::

    sys.constants = {lB: 0.2, # meters
                     h: 0.1, # meters
                     w: 0.2, # meters
                     mA: 0.01, # kilograms
                     mB: 0.1, # kilograms
                     g: 9.81} # meters per second squared

.. jupyter-execute::

    sys.initial_conditions = {theta: np.deg2rad(45),
                              phi: np.deg2rad(0.5),
                              omega: 0,
                              alpha: 0}

.. jupyter-execute::

    sys.times = np.linspace(0.0, 10.0, num=300)

The trajectories of the states are found with the ``integrate`` method.

.. jupyter-execute::

    x = sys.integrate()

The angles can be plotted to see how they change with respect to time
given the initial conditions.

.. jupyter-execute::

    def plot():
        plt.figure()
        plt.plot(sys.times, np.rad2deg(x[:, :2]))
        plt.legend([sm.latex(s, mode='inline') for s in sys.coordinates])

    plot()


Chaotic Behavior
================

Now change the intial condition of the plat angle just slighty to see if
the behvior of the system is similar.

.. jupyter-execute::

    sys.initial_conditions[phi] = np.deg2rad(1.0)
    x = sys.integrate()
    plot()

Seems all good, very similar behavior. But now set the rod angle to
:math:`90^\circ` and try the same slight change in plate angle.

.. jupyter-execute::

    sys.initial_conditions[theta] = np.deg2rad(90)
    sys.initial_conditions[phi] = np.deg2rad(0.5)
    x = sys.integrate()
    plot()

First note that the plate behaves wildly. What happens when the initial
plate angle is altered slightly.

.. jupyter-execute::

    sys.initial_conditions[phi] = np.deg2rad(1.0)
    x = sys.integrate()
    plot()

The behavior does not look similar to the previous simulation. This is
an example of chaotic behavior. The plate angle can not be reliably
predicted because slight changes in the initial conditions cause the
behavior of the system to vary widely.

Visualization
=============

Finally, the system can be animated by attached a cylinder and a plane
shape to the rigid bodies. To properly align the coordinate axes of the
shapes with the bodies, simple rotations are used.

.. jupyter-execute::

    rod_shape = Cylinder(2 * lA, 0.005, color='red', name='rod')
    plate_shape = Plane(w, h, color='blue', name='plate')

    v1 = VisualizationFrame('rod',
                            A.orientnew('rod', 'Axis', (sm.pi / 2, A.x)),
                            Ao,
                            rod_shape)

    v2 = VisualizationFrame('plate',
                            B.orientnew('plate', 'Body', (sm.pi / 2, sm.pi / 2, 0), 'XZX'),
                            Bo,
                            plate_shape)

    scene = Scene(N, No, v1, v2, system=sys)

The following method opens up a simple gui that shows a 3D animatoin of
the system.

.. jupyter-execute::

    scene.display_jupyter(axes_arrow_length=1.0)
