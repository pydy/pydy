#This one is compatible with SymPy 0.7.4.1

from sympy import symbols, sin, cos, tan
from sympy.physics.mechanics import *
x, y, q1, q2, q3, vx, vy, u1, u2, u3  = \
        dynamicsymbols('x y q1 q2 q3 xd yd u1 u2 u3')
xd, yd, q1d, q2d, q3d, vxd, vyd, u1d, u2d, u3d = \
        dynamicsymbols('x y q1 q2 q3 xd yd u1 u2 u3', 1)
r, m, g = symbols('r m g')
mechanics_printing()

N = ReferenceFrame('N')
Y = N.orientnew('Y', 'Axis', [q1, N.z])
L = Y.orientnew('L', 'Axis', [q2, Y.x])
R = L.orientnew('R', 'Axis', [q3, L.y])
w_R_N_qd = R.ang_vel_in(N)
R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

No = Point('No')
No.set_vel(N, 0)
Dmc = No.locatenew('Dmc', x * N.x + y * N.y + r * L.z)
Dmc.set_vel(N, xd * N.x + yd * N.y + L.ang_vel_in(N) ^ (r * L.z));
# Contact point.
C = Dmc.locatenew('C', -r * L.z)
C.v2pt_theory(Dmc, N, R)
#C = No.locatenew('C', x * N.x + y * N.y)
#C.set_vel(N, xd * N.x + yd * N.y)
#Dmc = C.locatenew('Dmc', r * L.z)
#Dmc.v2pt_theory(C, N, R)
I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
mprint(I)

kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L] + [vx - xd, vy - yd]

velocity_constraints = [C.vel(N) & Y.x , # slip in the direction of rolling.
                        C.vel(N) & Y.y, # perpendicular to the dir. of rolling.
                        ]

ForceList = [(Dmc, - m * g * Y.z)]
BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
BodyList = [BodyD]

KM = KanesMethod(N,
        q_ind=[x, y, q1, q2, q3],
        u_ind=[u1, u2, u3],
        kd_eqs=kd,
        u_dependent=[vx, vy],
        velocity_constraints=velocity_constraints)
(fr, frstar) = KM.kanes_equations(ForceList, BodyList)
MM = KM.mass_matrix
forcing = KM.forcing
rhs = MM.inv() * forcing
kdd = KM.kindiffdict()
rhs = rhs.subs(kdd)
rhs.simplify()
mprint(rhs)
