import sympy as sm
import sympy.physics.mechanics as me

q = me.dynamicsymbols('q:10')
u = me.dynamicsymbols('u:10')

pendulum_length = sm.symbols('l')
pendulum_base_mass = sm.symbols('m1')
pendulum_bob_mass = sm.symbols('m2')

kin_diffy_qs = [q_i.diff() - u_i for q_i, u_i in zip(q, u)]

newtonian = me.ReferenceFrame('N')
main_shaft = newtonian.orientnew('A', 'Quaternion', q[:4])
plane = main_shaft.orientnew('B', 'Axis', (q[5], main_shaft.z))
pendulum = plane.orientnew('C', 'Axis', (q[9], plane.y))

spherical_vel = main_shaft.ang_vel_in(newtonian)

main_shaft.set_ang_vel(newtonian, spherical_vel.subs({q_i.diff(): u_i for q_i, u_i in zip(q, u)}))
plane.set_ang_vel(main_shaft, u[5] * main_shaft.z)
pendulum.set_ang_vel(plane, u[9] * plane.y)

O = me.Point('O')
plane_mass_center = O.locatenew('P1', q[4] * main_shaft.z)
box_mass_center = plane_mass_center.locatenew('P2', q[6] * plane.x + q[7] *
                                              plane.z)
pendulum_base = plane_mass_center.locatenew('P3', q[8] * plane.x)
pendulum_bob = pendulum_base.locatenew('P4', l * pendulum.z)

O.set_vel(newtonian, 0)
plane_mass_center.set_vel(main_shaft, u[4] * main_shaft.z)
box_mass_center.set_vel(plane, u[6] * plane.x + u[7] * plane.z)
pendulum_base.set_vel(plane, u[8] * plane.x)
pendulum_bob.v2pt_theory(pendulum_base, newtonian, pendulum)

# forces

# gravitational
plane_mass * g * newtonian.z, p1
box_mass * g * newtonian.z, p2
pendulum_base_mass * g * newtonian.z, p3
pendulum_bob_mass * g * newtonian.z, p4

# spring forces
-k4 * q[4] * main_shaft.z, p1

-k6 * q[6] * plane.x - k7 * q[7] * plane.z, p2
k6 * q[6] * plane.x + k7 * q[7] * plane.z, p1

-k8 * q[8] * plane.x, p3
k8 * q[8] * plane.x, p1

# inertia

# TODO : Create basic inertia for a plane and block from standard formulas
# from dimensions and mass.

# bodies and particles

me.RigidBody('plane', p1, plane, plane_mass, plane_central_inertia)
me.RigidBody('block', p2, plane, box_mass, box_central_inertia)
me.Particle('pend_base', p3, pendulum_base_mass)
me.Particle('pend_bob', p4, pendulum_bob_mass)
