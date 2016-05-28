# This code requires sympy 1.0 to run

from sympy import *
from sympy.physics.mechanics import LagrangesMethod, Lagrangian
from sympy.physics.mechanics import ReferenceFrame, Particle, Point
from sympy.physics.mechanics import dynamicsymbols

# System state variables
theta = dynamicsymbols('theta')
thetad = dynamicsymbols('theta', 1)

# Other system variables
m, l, g = symbols('m l g')

# Set up the reference frames
# Reference frame A set up in the plane perpendicular to the page containing
# segment OP
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [theta, N.z])

# Set up the points and particles
O = Point('O')
P = O.locatenew('P', l * A.x)

Pa = Particle('Pa', P, m)

# Set up velocities
A.set_ang_vel(N, thetad * N.z)
O.set_vel(N, 0)
P.v2pt_theory(O, N, A)

# Set up the lagrangian
L = Lagrangian(N, Pa)

# Create the list of forces acting on the system
fl = [(P, g * m * N.x)]

# Create the equations of motion using lagranges method
l = LagrangesMethod(L, [theta], forcelist=fl, frame=N)

pprint(l.form_lagranges_equations())
