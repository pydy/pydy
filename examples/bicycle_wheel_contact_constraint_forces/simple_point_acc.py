# This scripts calculates the body fixed acceleration of the IMU attached on
# the rear frame ignoring any pitch degree of freedom.
import sympy as sm
import sympy.physics.mechanics as me

# bx, by, bz are the measures from C (contact point) to P (imu point on rear
# frame) along the B basis vectors.
g, bx, by, bz = sm.symbols('g, bx, by, bz')
q1, q2, q3, q4 = me.dynamicsymbols('q1, q2, q3, q4')
u1, u2, u3, u4 = me.dynamicsymbols('u1, u2, u3, u4')
u1p, u2p, u3p, u4p = me.dynamicsymbols('u1p, u2p, u3p, u4p')

N = me.ReferenceFrame('N')  # inertial frame
A = me.ReferenceFrame('A')  # yaw frame
B = me.ReferenceFrame('B')  # roll frame

# yaw then roll
A.orient(N, 'Axis', (q3, N.z))
B.orient(A, 'Axis', (q4, A.x))

A.set_ang_vel(N, u3*N.z)
B.set_ang_vel(A, u4*A.x)

# origin
O = me.Point('O')
O.set_vel(N, 0)

# contact point
C = O.locatenew('C', q1*N.x + q2*N.y)
C.set_vel(N, u1*N.x + u2*N.y)
C.set_acc(N, u1p*N.x + u2p*N.y)

# imu point
P = C.locatenew('P', bx*B.x + by*B.y + bz*B.z)
P.v2pt_theory(C, N, B)
P.a2pt_theory(C, N, B)

deriv_rpl = {
    u1.diff(): u1p,
    u2.diff(): u2p,
    u3.diff(): u3p,
    u4.diff(): u4p,
    q1.diff(): u1,
    q2.diff(): u2,
    q3.diff(): u3,
    q4.diff(): u4,
}

# body fixed acceleration in N that the IMU outputs (does not gravity)
N_A_P = P.acc(N).to_matrix(B).xreplace(deriv_rpl)

B_P = sm.trigsimp(N_A_P.jacobian([u4p, u3p]))
Y_accP = sm.trigsimp(-N_A_P.xreplace({u4p: 0, u3p: 0}))
