#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.3 from Kane 1985."""

from __future__ import division
from sympy import Matrix, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols


R_SQ, R_QE, R_EM, R_Ee, R_Mm = symbols('R_SQ R_QE R_EM R_Ee R_Mm', positive=True)
omega_E, omega_M, omega_e, omega_m = symbols('ω_E ω_M ω_e ω_m', positive=True)

symbol_values = {R_SQ: 4.5e5, R_QE: 1.5e11, R_EM: 4.0e8,
                 R_Ee: 7.0e6, R_Mm: 2.0e6,
                 omega_E: 2e-7, omega_M: 24e-7,
                 omega_e: 12e-4, omega_m:10e-4}

# reference frames
S = ReferenceFrame('S')
Q = S.orientnew('Q', 'axis', [0, S.x])
E = Q.orientnew('E', 'axis', [0, S.x])
M = E.orientnew('M', 'axis', [0, S.x])
frames = [S, Q, E, M]

pS = Point('S')
pS.set_acc(S, 0)
pQ = Point('Q')
pQ.set_acc(Q, 0)
pE = Point('E')
pE.set_acc(E, 0)
pM = Point('M')
pM.set_acc(M, 0)
pe = Point('e')
pm = Point('m')
points = [pS, pQ, pE, pM, pe, pm]

# v = ω*R, a = ω**2 * R
pQ.set_acc(S, omega_E**2 * R_SQ * S.x)
pE.set_acc(Q, omega_E**2 * R_QE * S.x)
pM.set_acc(E, omega_M**2 * R_EM * S.x)
pe.set_acc(E, omega_e**2 * R_Ee * S.x)
pm.set_acc(M, omega_m**2 * R_Mm * S.x)

# v_p_A = v_p_B + cross(ω_B_A, r)
#       = v_p_B, since angular vel is zero
# a_p_A = a_p_B + cross(ω_B_A, v_p_A)
#       = a_p_B
# and a_p_A = a_p_Bbar + a_p_B + 2*cross(ω_B_A, v_p_B)
#           = a_p_Bbar + a_p_B
pm.set_acc(E, pm.acc(M) + pM.acc(E))
pm.set_acc(Q, pm.acc(E) + pE.acc(Q))
pm.set_acc(S, pm.acc(Q) + pQ.acc(S))

pe.set_acc(M, pe.acc(E) + pM.acc(E))
pe.set_acc(Q, pe.acc(E) + pE.acc(Q))
pe.set_acc(S, pe.acc(Q) + pQ.acc(S))

pM.set_acc(Q, pM.acc(E) + pE.acc(Q))
pM.set_acc(S, pM.acc(Q) + pQ.acc(S))

pE.set_acc(M, pE.acc(E) + pM.acc(E))
pE.set_acc(S, pE.acc(Q) + pQ.acc(S))

pQ.set_acc(E, pQ.acc(Q) + pE.acc(Q))
pQ.set_acc(M, pQ.acc(E) + pM.acc(E))

pS.set_acc(Q, pS.acc(S) + pQ.acc(S))
pS.set_acc(E, pS.acc(Q) + pE.acc(Q))
pS.set_acc(M, pS.acc(E) + pM.acc(E))

#print('acc in frame\t{0}\t{1}\t{2}\t{3}'.format(*frames))
#for p in points:
#    print('point {0}:\t{1:0.3g}\t{2:0.3g}\t{3:0.3g}\t{4:0.3g}'.format(
#    #print('point {0}:\t{1}\t{2}\t{3}\t{4}'.format(
#            p, *map(lambda x: float(p.acc(x).magnitude().subs(symbol_values)),
#                    frames)))


idx_Q = frames.index(Q)
print('acc in frame ratio\t{0}\t{1}\t{2}'.format(S, E, M))
acc_ratios = Matrix.zeros(4, 3)
for i, p in enumerate(points[2:]):
    acc_values = map(lambda x: p.acc(x).magnitude().subs(symbol_values),
                     frames)
    a_p_Q = acc_values[idx_Q]
    acc_values = [float(x / a_p_Q)
                  for x in acc_values[:idx_Q] + acc_values[idx_Q + 1:]]
    acc_ratios[i, :] = Matrix(acc_values).T
    print('object {0}:\t\t{1:0.3g}\t{2:0.3g}\t{3:0.3g}'.format(
            p, *acc_values))
print('Approximately Newtonian reference frames have a near 1 ratio.')

min_ratio = 0.9 # minimum value if frame is approximately Newtonian.
assert acc_ratios[0, 0] >= min_ratio
assert acc_ratios[0, 1] < min_ratio
assert acc_ratios[0, 2] < min_ratio
assert acc_ratios[1, 0] >= min_ratio
assert acc_ratios[1, 1] < min_ratio
assert acc_ratios[1, 2] < min_ratio
assert acc_ratios[2, 0] >= min_ratio
assert acc_ratios[2, 1] >= min_ratio
assert acc_ratios[2, 2] >= min_ratio
assert acc_ratios[3, 0] >= min_ratio
assert acc_ratios[3, 1] >= min_ratio
assert acc_ratios[3, 2] >= min_ratio
