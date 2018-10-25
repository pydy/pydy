"""
Example to solve Kinematics of 2DoF Scara Arm.
"""
# Funções das Bibliotecas Utilizadas
from sympy import symbols, trigsimp, pprint
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector import ReferenceFrame, Vector
from sympy.physics.vector import time_derivative

# Variáveis Simbólicas
THETA_1, THETA_2 = dynamicsymbols('theta_1 theta_2')
L_1, L_2 = symbols('l_1 l_2', positive=True)

# Referenciais
# Referencial Parado
B0 = ReferenceFrame('B0')
B1 = ReferenceFrame('B1')
# Referencial móvel: theta_1 em relação a B0.z
B1.orient(B0, 'Axis', [THETA_1, B0.z])
B2 = ReferenceFrame('B2')
# Referencial móvel: theta_2 em relação a B1.z
B2.orient(B1, 'Axis', [THETA_2, B1.z])

# Vetores Posição entre os Pontos
# Vetor Nulo
B0_R_OA = Vector(0)
# Vetor que liga os pontos A e B expresso no referencial móvel B1
B1_R_AB = L_1 * B1.x
# Vetor que liga os pontos B e C expresso no referencial móvel B2
B2_R_BC = L_2 * B2.x

# Cinemática do ponto A em relação ao referencial B0
R_A = B0_R_OA
V_A = time_derivative(R_A, B0)
A_A = time_derivative(V_A, B0)

# Cinemática do ponto B em relação ao referencial B0
R_B = R_A + B1_R_AB.express(B0)
V_B = time_derivative(R_B, B0)
A_B = time_derivative(V_B, B0)

# Cinemática do ponto C em relação ao referencial B0
R_C = R_B.express(B0) + B2_R_BC.express(B0)
V_C = (time_derivative(R_C, B0))
A_C = (time_derivative(V_C, B0))

# Simplificação dos Resultados
R_A = (R_A.to_matrix(B0)).applyfunc(trigsimp)
V_A = (V_A.to_matrix(B0)).applyfunc(trigsimp)
A_A = (A_A.to_matrix(B0)).applyfunc(trigsimp)
R_B = (R_B.to_matrix(B0)).applyfunc(trigsimp)
V_B = (V_B.to_matrix(B0)).applyfunc(trigsimp)
A_B = (A_B.to_matrix(B0)).applyfunc(trigsimp)
R_C = (R_C.to_matrix(B0)).applyfunc(trigsimp)
V_C = (V_C.to_matrix(B0)).applyfunc(trigsimp)
A_C = (A_C.to_matrix(B0)).applyfunc(trigsimp)

# Resultados de C
pprint(R_C)
pprint(V_C)
pprint(A_C)
