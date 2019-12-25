===========================
Three Link Conical Pendulum
===========================

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`three-link-conical-pendulum` or Jupyter notebook:
   :jupyter-download:notebook:`three-link-conical-pendulum`.

This example shows how to simulate a three link conical compound pendulum made
up of cylindrical bars connecting particles at the joints. The is an example of
a 3D holonomic system. Wikipedia provides a basic description of a `conical
pendulum`_.

.. _conical pendulum: https://en.wikipedia.org/wiki/Conical_pendulum

Derive the Equations of Motion
==============================

.. jupyter-execute::

   from sympy import symbols
   import sympy.physics.mechanics as me

   print("Defining the problem.")

   # The conical pendulum will have three links and three bobs.
   n = 3

   # Each link's orientation is described by two spaced fixed angles: alpha and
   # beta.

   # Generalized coordinates
   alpha = me.dynamicsymbols('alpha:{}'.format(n))
   beta = me.dynamicsymbols('beta:{}'.format(n))

   # Generalized speeds
   omega = me.dynamicsymbols('omega:{}'.format(n))
   delta = me.dynamicsymbols('delta:{}'.format(n))

   # At each joint there are point masses (i.e. the bobs).
   m_bob = symbols('m:{}'.format(n))

   # Each link is modeled as a cylinder so it will have a length, mass, and a
   # symmetric inertia tensor.
   l = symbols('l:{}'.format(n))
   m_link = symbols('M:{}'.format(n))
   Ixx = symbols('Ixx:{}'.format(n))
   Iyy = symbols('Iyy:{}'.format(n))
   Izz = symbols('Izz:{}'.format(n))

   # Acceleration due to gravity will be used when prescribing the forces
   # acting on the links and bobs.
   g = symbols('g')

   # Now defining an inertial reference frame for the system to live in. The Y
   # axis of the frame will be aligned with, but opposite to, the gravity
   # vector.

   I = me.ReferenceFrame('I')

   # Three reference frames will track the orientation of the three links.

   A = me.ReferenceFrame('A')
   A.orient(I, 'Space', [alpha[0], beta[0], 0], 'ZXY')

   B = me.ReferenceFrame('B')
   B.orient(A, 'Space', [alpha[1], beta[1], 0], 'ZXY')

   C = me.ReferenceFrame('C')
   C.orient(B, 'Space', [alpha[2], beta[2], 0], 'ZXY')

   # Define the kinematical differential equations such that the generalized
   # speeds equal the time derivative of the generalized coordinates.
   kinematic_differentials = []
   for i in range(n):
      kinematic_differentials.append(omega[i] - alpha[i].diff())
      kinematic_differentials.append(delta[i] - beta[i].diff())

   # The angular velocities of the three frames can then be set.
   A.set_ang_vel(I, omega[0] * I.z + delta[0] * I.x)
   B.set_ang_vel(I, omega[1] * I.z + delta[1] * I.x)
   C.set_ang_vel(I, omega[2] * I.z + delta[2] * I.x)

   # The base of the pendulum will be located at a point O which is stationary
   # in the inertial reference frame.
   O = me.Point('O')
   O.set_vel(I, 0)

   # The location of the bobs (at the joints between the links) are created by
   # specifiying the vectors between the points.
   P1 = O.locatenew('P1', -l[0] * A.y)
   P2 = P1.locatenew('P2', -l[1] * B.y)
   P3 = P2.locatenew('P3', -l[2] * C.y)

   # The velocities of the points can be computed by taking advantage that
   # pairs of points are fixed on the referene frames.
   P1.v2pt_theory(O, I, A)
   P2.v2pt_theory(P1, I, B)
   P3.v2pt_theory(P2, I, C)
   points = [P1, P2, P3]

   # Now create a particle to represent each bob.
   Pa1 = me.Particle('Pa1', points[0], m_bob[0])
   Pa2 = me.Particle('Pa2', points[1], m_bob[1])
   Pa3 = me.Particle('Pa3', points[2], m_bob[2])
   particles = [Pa1, Pa2, Pa3]

   # The mass centers of each link need to be specified and, assuming a
   # constant density cylinder, it is equidistance from each joint.
   P_link1 = O.locatenew('P_link1', -l[0] / 2 * A.y)
   P_link2 = P1.locatenew('P_link2', -l[1] / 2 * B.y)
   P_link3 = P2.locatenew('P_link3', -l[2] / 2 * C.y)

   # The linear velocities can be specified the same way as the bob points.
   P_link1.v2pt_theory(O, I, A)
   P_link2.v2pt_theory(P1, I, B)
   P_link3.v2pt_theory(P2, I, C)

   points_rigid_body = [P_link1, P_link2, P_link3]

   # The inertia tensors for the links are defined with respect to the mass
   # center of the link and the link's reference frame.
   inertia_link1 = (me.inertia(A, Ixx[0], Iyy[0], Izz[0]), P_link1)
   inertia_link2 = (me.inertia(B, Ixx[1], Iyy[1], Izz[1]), P_link2)
   inertia_link3 = (me.inertia(C, Ixx[2], Iyy[2], Izz[2]), P_link3)

   # Now rigid bodies can be created for each link.
   link1 = me.RigidBody('link1', P_link1, A, m_link[0], inertia_link1)
   link2 = me.RigidBody('link2', P_link2, B, m_link[1], inertia_link2)
   link3 = me.RigidBody('link3', P_link3, C, m_link[2], inertia_link3)
   links = [link1, link2, link3]

   # The only contributing forces to the system is the force due to gravity
   # acting on each particle and body.
   forces = []

   for particle in particles:
      mass = particle.mass
      point = particle.point
      forces.append((point, -mass * g * I.y))

   for link in links:
      mass = link.mass
      point = link.masscenter
      forces.append((point, -mass * g * I.y))

   # Make a list of all the particles and bodies in the system.
   total_system = links + particles

   # Lists of all generalized coordinates and speeds.
   q = alpha + beta
   u = omega + delta

   # Now the equations of motion of the system can be formed.
   print("Generating equations of motion.")
   kane = me.KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kinematic_differentials)
   fr, frstar = kane.kanes_equations(total_system, loads=forces)
   print("Derivation complete.")

Simulate the System
===================

.. jupyter-execute::

   # external
   from numpy import radians, linspace, hstack, zeros, ones
   from scipy.integrate import odeint
   from pydy.codegen.ode_function_generators import generate_ode_function

   param_syms = []
   for par_seq in [l, m_bob, m_link, Ixx, Iyy, Izz, (g,)]:
      param_syms += list(par_seq)

   # All of the links and bobs will have the same numerical values for the
   # parameters.

   link_length = 10.0  # meters
   link_mass = 10.0  # kg
   link_radius = 0.5  # meters
   link_ixx = 1.0 / 12.0 * link_mass * (3.0 * link_radius**2 + link_length**2)
   link_iyy = link_mass * link_radius**2
   link_izz = link_ixx

   particle_mass = 5.0  # kg
   particle_radius = 1.0  # meters

   # Create a list of the numerical values which have the same order as the
   # list of symbolic parameters.
   param_vals = [link_length for x in l] + \
               [particle_mass for x in m_bob] + \
               [link_mass for x in m_link] + \
               [link_ixx for x in list(Ixx)] + \
               [link_iyy for x in list(Iyy)] + \
               [link_izz for x in list(Izz)] + \
               [9.8]

   # A function that evaluates the right hand side of the set of first order
   # ODEs can be generated.
   print("Generating numeric right hand side.")
   right_hand_side = generate_ode_function(kane.forcing_full, q, u, param_syms,
                                          mass_matrix=kane.mass_matrix_full,
                                          generator='cython')

   # To simulate the system, a time vector and initial conditions for the
   # system's states is required.
   duration = 10.0
   fps = 60.0
   t = linspace(0.0, duration, num=int(duration*fps))
   x0 = hstack((ones(6) * radians(10.0), zeros(6)))

   print("Integrating equations of motion.")
   state_trajectories = odeint(right_hand_side, x0, t, args=(dict(zip(param_syms,
                                                                     param_vals)),
                                                            ))
   print("Integration done.")

Visualize the System
====================

.. jupyter-execute::

   # external
   from pydy.viz.shapes import Cylinder, Sphere
   from pydy.viz.scene import Scene
   from pydy.viz.visualization_frame import VisualizationFrame

   # A cylinder will be attached to each link and a sphere to each bob for the
   # visualization.

   viz_frames = []

   for i, (link, particle) in enumerate(zip(links, particles)):

      link_shape = Cylinder(name='cylinder{}'.format(i),
                           radius=link_radius,
                           length=link_length,
                           color='red')

      viz_frames.append(VisualizationFrame('link_frame{}'.format(i), link,
                                          link_shape))

      particle_shape = Sphere(name='sphere{}'.format(i),
                              radius=particle_radius,
                              color='blue')

      viz_frames.append(VisualizationFrame('particle_frame{}'.format(i),
                                          link.frame,
                                          particle,
                                          particle_shape))

   # Now the visualization frames can be passed in to create a scene.
   scene = Scene(I, O, *viz_frames)

   # Provide the data to compute the trajectories of the visualization frames.
   scene.times = t
   scene.constants = dict(zip(param_syms, param_vals))
   scene.states_symbols = q + u
   scene.states_trajectories = state_trajectories

   scene.display_jupyter()
