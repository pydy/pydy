import sympy as sp
import sympy.physics.mechanics as me

# declare the coordinates and speeds and their derivatives #
# theta : angle of the rod
# phi : angle of the plate
#!/usr/bin/env python

# omega : angular speed of the rod
# alpha : angular speed of the plate
theta, phi, omega, alpha = me.dynamicsymbols('theta phi omega alpha')
thetad, phid, omegad, alphad = me.dynamicsymbols('theta phi omega alpha', 1)

# declare the constants #
# gravity
gravity = sp.symbols('g')
# center of mass length, mass and  moment of inertia of the slender rod
lA, mA, IA = sp.symbols('lA mA IA')
# center of mass length, mass and moment of inertia of the plate
lB, mB, IB1, IB2, IB3 = sp.symbols('lB mB IB1 IB2 IB3')

# create a Newtonian reference frame
N = me.ReferenceFrame('N', indices=['1', '2', '3'])

# orientations #
# the rod rotates with respect to the Newtonian reference frame about the 2
# axis
A = N.orientnew('A', 'Axis', [theta, N['2']])
# the plate rotates about the rod's primay axis
B = A.orientnew('B', 'Axis', [phi, N['3']])

# positions #
# original of the Newtonian reference frame
No = me.Point('No')
# create a point for the mass centers of the two bodies
Ao = me.Point('Ao')
Bo = me.Point('Bo')
# define the positions of the mass centers relative to the Newtonian origin
Ao.set_pos(No, lA * A['3'])
Bo.set_pos(No, lB * A['3'])

# angular velocities and accelerations
A.set_ang_vel(N, omega * N['2'])
B.set_ang_vel(A, alpha * A['3'])

A.set_ang_acc(N, A.ang_vel_in(N).dt(N))
B.set_ang_acc(A, B.ang_vel_in(A).dt(A))

# linear velocities and accelerations
N.set_vel(N, 0)
B.set_vel(N, omega * l * a.x)
B.a2pt(N, n, b)

I = inertia(b, I1, I2, I3)

kd = [omega-thetad, alpha-phid]
forcelist = [(B, -n.y*gravity*m)]
plate = RigidBody()
plate.mc = B
plate.inertia = (I, B)
plate.frame = B
print b.dcm(n)
