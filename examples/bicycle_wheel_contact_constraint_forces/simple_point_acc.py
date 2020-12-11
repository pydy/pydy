# This scripts calculates the body fixed acceleration of the two IMUs attached
# to the rear and front frames of the bicycle, ignoring any pitch degree of
# freedom.
import sympy as sm
import sympy.physics.mechanics as me

lam = sm.symbols('lambda')  # steer axis tilt
# px, py, pz are the measures from C (contact point) to P (imu point on rear
# frame) along the B basis vectors.
px, py, pz = sm.symbols('p_x, p_y, p_z')
# hx, hy, hz are the measures from C (contact point) to H (steer axis point)
# along the B basis vectors.
hx, hy, hz = sm.symbols('h_x, h_y, h_z')
# qx, qy, qz are the measures from H (steer axis point) to Q (imu point on
# front frame) along the D basis vectors.
qx, qy, qz = sm.symbols('q_x, q_y, q_z')

# generalized coordinates, speeds, and accelerations
q1, q2, q3, q4, q7 = me.dynamicsymbols('q1, q2, q3, q4, q7')
u1, u2, u3, u4, u7 = me.dynamicsymbols('u1, u2, u3, u4, u7')
u1p, u2p, u3p, u4p, u7p = me.dynamicsymbols('u_{1p}, u_{2p}, u_{3p}, u_{4p}, u_{7p}')

# reference frames
N = me.ReferenceFrame('N')  # inertial frame
A = me.ReferenceFrame('A')  # yaw frame
B = me.ReferenceFrame('B')  # roll frame
D = me.ReferenceFrame('D')  # steer frame

# yaw then roll
A.orient(N, 'Axis', (q3, N.z))
B.orient(A, 'Axis', (q4, A.x))
# pitch through constant lambda and steer
D.orient(B, 'Body', (lam, q7, 0), rot_order='yzy')

# angular velocity of simple rotations
A.set_ang_vel(N, u3*N.z)
B.set_ang_vel(A, u4*A.x)
D.set_ang_vel(A, u7*D.z)

# origin
O = me.Point('O')
O.set_vel(N, 0)

# contact point
C = O.locatenew('C', q1*N.x + q2*N.y)
C.set_vel(N, u1*N.x + u2*N.y)
C.set_acc(N, u1p*N.x + u2p*N.y)

# rear imu point
P = C.locatenew('P', px*B.x + py*B.y + pz*B.z)
P.v2pt_theory(C, N, B)
P.a2pt_theory(C, N, B)

# steer axis point
H = C.locatenew('H', hx*B.x + hy*B.y + hz*B.z)
H.v2pt_theory(C, N, B)
H.a2pt_theory(C, N, B)

# front imu point
Q = H.locatenew('Q', qx*D.x + qy*D.y + qz*D.z)
Q.v2pt_theory(H, N, D)
Q.a2pt_theory(H, N, D)

deriv_rpl = {
    u1.diff(): u1p,
    u2.diff(): u2p,
    u3.diff(): u3p,
    u4.diff(): u4p,
    u7.diff(): u7p,
    q1.diff(): u1,
    q2.diff(): u2,
    q3.diff(): u3,
    q4.diff(): u4,
    q7.diff(): u7,
}

# body fixed acceleration in N that the P IMU outputs (does not include gravity
# like the IMU does)
N_A_P = P.acc(N).to_matrix(B).xreplace(deriv_rpl)
B_P = N_A_P.jacobian([u4p, u3p])
Y_accP = -N_A_P.xreplace({u4p: 0, u3p: 0})

# body fixed acceleration in N that the Q IMU outputs (does not include gravity
# like the IMU does)
N_A_Q = Q.acc(N).to_matrix(D).xreplace(deriv_rpl)
B_Q = N_A_Q.jacobian([u4p, u3p, u7p])
Y_accQ = -N_A_Q.xreplace({u4p: 0, u3p: 0, u7p: 0})
