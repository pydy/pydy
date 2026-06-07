"""Example demonstrating how to obtain the equations of motion of a simple
pendulum with Lagrange's method."""

from sympy import symbols, pprint
from sympy.physics.mechanics import LagrangesMethod, Lagrangian
from sympy.physics.mechanics import ReferenceFrame, Particle, Point
from sympy.physics.mechanics import dynamicsymbols

# The pendulum angle is the single generalized coordinate.
theta = dynamicsymbols('theta')

# The bob has a mass m, the pendulum lenth is l, and acceleration due to
# gravity is g.
m, l, g = symbols('m, l, g')

# Rotate the reference frame A relative to the inertial reference frame N.
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', (theta, N.z))

# Locate the bob relative to the hinge, which is fixed in N.
O = Point('O')
O.set_vel(N, 0)
P = O.locatenew('P', l*A.x)

# The velocity of the bob can be calculated.
pprint(P.vel(N))

# Create the Lagrangian for the system by ensurnig the particle has kinetic and
# potential energy.
Pa = Particle('Pa', P, m)
Pa.potential_energy = -m*g*P.pos_from(O).dot(N.x)
L = Lagrangian(N, Pa)

# Create the equations of motion using Lagrange's method.
lm = LagrangesMethod(L, [theta], frame=N)
pprint(lm.form_lagranges_equations())
