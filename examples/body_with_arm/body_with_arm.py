from sympy import symbols
from sympy.physics.mechanics import *

q = dynamicsymbols('q:6')
qd = dynamicsymbols('q:6', 1)
u = dynamicsymbols('u:6')
ud = dynamicsymbols('u:6', 1)
x, theta = dynamicsymbols('x theta')
xd, thetad = dynamicsymbols('x theta', 1)
l = symbols('l')


N = ReferenceFrame('N')
A = N.orientnew('A', 'Body', q[0:3], 313)
A.set_ang_vel(N, u[0] * A.x + u[1] * A.y + u[2] * A.z)
B = A.orientnew('B', 'Axis', (theta, A.x))

O = Point('O')
O.set_vel(N, u[3] * A.x + u[4] * A.y + u[5] * A.z)
P = O.locatenew('P', -l * B.z + x * A.y)
P.set_vel(B, x * B.y)

P.v1pt_theory(O, N, B)
P.a1pt_theory(O, N, B)
