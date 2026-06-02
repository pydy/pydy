"""
Example to solve Dynamics of 2DoF Scara Arm.
"""
# Funções e Bibliotecas Utilizadas
from sympy import symbols, pprint
from sympy.physics.mechanics import (
    dynamicsymbols,
    ReferenceFrame,
    Point,
    RigidBody,
    init_vprinting,
    Lagrangian,
    LagrangesMethod)
from sympy.physics.mechanics.functions import inertia
init_vprinting()

# Variáveis Simbólicas
THETA_1, THETA_2 = dynamicsymbols('theta_1 theta_2')
DTHETA_1, DTHETA_2 = dynamicsymbols('theta_1 theta_2', 1)
TAU_1, TAU_2 = symbols('tau_1 tau_2')
L_1, L_2 = symbols('l_1 l_2', positive=True)
R_1, R_2 = symbols('r_1 r_2', positive=True)
M_1, M_2, G = symbols('m_1 m_2 g')
I_1_ZZ, I_2_ZZ = symbols('I_{1zz}, I_{2zz}')

# Referenciais
# Referencial Inercial
B0 = ReferenceFrame('B0')
# Referencial móvel: theta_1 em relação a B0.z
B1 = ReferenceFrame('B1')
B1.orient(B0, 'Axis', [THETA_1, B0.z])
# Referencial móvel: theta_2 em relação a B1.z
B2 = ReferenceFrame('B2')
B2.orient(B1, 'Axis', [THETA_2, B1.z])

# Pontos e Centros de Massa
O = Point('O')
O.set_vel(B0, 0)
A = Point('A')
A.set_pos(O, L_1 * B1.x)
A.v2pt_theory(O, B0, B1)
CM_1 = Point('CM_1')
CM_1.set_pos(O, R_1 * B1.x)
CM_1.v2pt_theory(O, B0, B1)
CM_2 = Point('CM_2')
CM_2.set_pos(A, R_2 * B2.x)
CM_2.v2pt_theory(O, B0, B2)

# Corpos Rígidos
I_1 = inertia(B1, 0, 0, I_1_ZZ)
# Elo 1
E_1 = RigidBody('Elo_1', CM_1, B1, M_1, (I_1, CM_1))
I_2 = inertia(B1, 0, 0, I_1_ZZ)
# Elo 2
E_2 = RigidBody('Elo_2', CM_2, B2, M_2, (I_2, CM_2))

# Energia Potencial
P_1 = -M_1 * G * B0.z
R_1_CM = CM_1.pos_from(O).express(B0)
E_1.potential_energy = R_1_CM.dot(P_1)
P_2 = -M_2 * G * B0.z
R_2_CM = CM_2.pos_from(O).express(B0).simplify()
E_2.potential_energy = R_2_CM.dot(P_2)

# Forças/Momentos Generalizados
FL = [(B1, TAU_1 * B0.z), (B2, TAU_2 * B0.z)]
# Método de Lagrange
L = Lagrangian(B0, E_1, E_2)
L = L.simplify()
pprint(L)
LM = LagrangesMethod(L, [THETA_1, THETA_2], frame=B0, forcelist=FL)

# Equações do Movimento
L_EQ = LM.form_lagranges_equations()
pprint(L_EQ)
# Equações Prontas para Solução Numérica
RHS = LM.rhs()
