import sympy as sm
import sympy.physics.mechanics as me

q1, q2 = me.dynamicsymbols('q1 q2')
q1d, q2d = me.dynamicsymbols('q1 q2', 1)
u1, u2 = me.dynamicsymbols('u1 u2')
u1d, u2d = me.dynamicsymbols('u1 u2', 1)
l, m, g = sm.symbols('l m g')

N = me.ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [q1, N.z])
B = N.orientnew('B', 'Axis', [q2, N.z])

A.set_ang_vel(N, u1 * N.z)
B.set_ang_vel(N, u2 * N.z)

O = me.Point('O')
P = O.locatenew('P', l * A.x)
R = P.locatenew('R', l * B.x)

O.set_vel(N, 0)
P.v2pt_theory(O, N, A)
R.v2pt_theory(P, N, B)

ParP = me.Particle('ParP', P, m)
ParR = me.Particle('ParR', R, m)

kd = [q1d - u1, q2d - u2]
FL = [(P, m * g * N.x), (R, m * g * N.x)]
BL = [ParP, ParR]

KM = me.KanesMethod(N, q_ind=[q1, q2], u_ind=[u1, u2], kd_eqs=kd)
KM.kanes_equations(BL, loads=FL)

kdd = KM.kindiffdict()
mass_matrix = KM.mass_matrix_full
forcing_vector = KM.forcing_full
qudots = mass_matrix.inv() * forcing_vector
qudots = qudots.subs(kdd)
