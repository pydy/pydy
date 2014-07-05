#This one is compatible with SymPy 0.7.4.1

from sympy import symbols, solve, Matrix, eye, zeros, sin
import sympy.physics.mechanics as me

"""
x : x coordinate of the ground contact path point in the ground plane
y : y coordinate of the ground contact path point in the ground plane

q1 : yaw
q2 : roll
q3 : pitch

"""

x, y, q1, q2, q3 = me.dynamicsymbols('x y q1 q2 q3')
vx, vy, u1, u2, u3 = me.dynamicsymbols('vx vy u1 u2 u3')

xd, yd, q1d, q2d, q3d = me.dynamicsymbols('x y q1 q2 q3', 1)
vxd, vyd, u1d, u2d, u3d = me.dynamicsymbols('vx vy u1 u2 u3', 1)

r, m, g = symbols('r m g')

me.mechanics_printing(pretty_print=True)

N = me.ReferenceFrame('N')
Y = N.orientnew('Y', 'Axis', [q1, N.z])  # Yaw
L = Y.orientnew('L', 'Axis', [q2, Y.x])  # Roll/Lean
R = L.orientnew('R', 'Axis', [q3, L.y])  # Pitch

# Specify the kinematic differential equations
w_R_N_qd = R.ang_vel_in(N)
R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
kd = [me.dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L] + [vx - xd, vy - yd]
qdots = solve(kd, (xd, yd, q1d, q2d, q3d))

# Global origin
No = me.Point('No')
No.set_vel(N, 0)

# Disc mass center
Ro = No.locatenew('Ro', x * N.x + y * N.y + r * L.z)
# L.ang_vel(N) will have qdots in the expression, so they will need to be
# subbed out with u's.
Ro.set_vel(N, (vx * N.x + vy * N.y + me.cross(L.ang_vel_in(N), r * L.z)).subs(qdots))

# Contact point.
C = Ro.locatenew('C', -r * L.z)
C.v2pt_theory(Ro, N, R)

print('Kinematical Differential Equations')
me.mprint(kd)

velocity_constraints = [
                        me.dot(C.vel(N), Y.x), # slip in the direction of rolling.
                        me.dot(C.vel(N),  Y.y), # perpendicular to the dir. of rolling.
                       ]

ForceList = [(Ro, -m * g * Y.z)]

I = me.inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
BodyR = me.RigidBody('BodyR', Ro, R, m, (I, Ro))

BodyList = [BodyR]

# We can solve the velocity constraints such that they are completely in
# terms of the independent u's.
vs = solve(velocity_constraints, (vx, vy))

# Now differentiate these expressions and substitute the qdots, this gives
# expressions for the derivatives of the dependent speeds in terms of the
# independent speeds.

#vdots = {k.diff(): v.diff().subs(qdots).simplify() for k, v in vs.items()}
vdots = {k.diff(): v.diff(me.dynamicsymbols._t).subs(qdots).simplify() for k, v in vs.items()}

KM = me.KanesMethod(N,
                    q_ind=[q1, q2, q3, x, y],
                    u_ind=[u1, u2, u3],
                    kd_eqs=kd,
                    u_dependent=[vx, vy],
                    velocity_constraints=velocity_constraints)

(fr, frstar) = KM.kanes_equations(ForceList, BodyList)

# So at this point it seems like the dependent speeds are not eliminated
# from the equations. Neither the dynamic or kinematic equations.

zero = fr + frstar

zero = zero.subs(vdots)
zero.simplify()

coordinates = [q1, q2, q3, x, y]
speeds = [u1, u2, u3]

F_kin = Matrix([qdots[c.diff()].subs(vs) for c in coordinates])
M_kin = eye(len(coordinates))

M_dyn = []
F_dyn = []
for j, s in enumerate(speeds):
    M_dyn_row = []
    reduced = zero[j].copy()
    for state in coordinates + speeds:
        coeff = zero[j].expand().coeff(state.diff())
        M_dyn_row.append(-coeff)
        reduced -= coeff * state.diff()
    F_dyn.append(reduced)
    M_dyn.append(M_dyn_row)

M_dyn = Matrix(M_dyn)
F_dyn = Matrix(F_dyn)

F_dyn.simplify()

M_kin = M_kin.row_join(zeros(len(coordinates), len(speeds)))

M = M_kin.col_join(M_dyn)
F = F_kin.col_join(F_dyn)

# F may not evaluate exactly as F in the method that derives without the
# ignorable coordinates.
