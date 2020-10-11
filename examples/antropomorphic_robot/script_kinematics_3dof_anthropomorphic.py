"""
Example for Anthropomorphic Arm.
"""
# Funções das Bibliotecas Utilizadas
from sympy import symbols, trigsimp, pprint
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector import ReferenceFrame, Vector
from sympy.physics.vector import time_derivative

# Variáveis Simbólicas
THETA_1, THETA_2, THETA_3 = dynamicsymbols('THETA_1 THETA_2 THETA_3')
L_1, L_2 = symbols('L_1 L_2', positive=True)

# Referenciais
# Referencial Parado
B0 = ReferenceFrame('B0')
# Referencial móvel: THETA_1 em relação a B0.y
B1 = ReferenceFrame('B1')
B1.orient(B0, 'Axis', [THETA_1, B0.y])
# Referencial móvel: THETA_2 em relação a B1.z
B2 = ReferenceFrame('B2')
B2.orient(B1, 'Axis', [THETA_2, B1.z])
# Referencial móvel: THETA_3 em relação a B2.z
B3 = ReferenceFrame('B3')
B3.orient(B2, 'Axis', [THETA_3, B2.z])

# Vetores Posição entre os Pontos
# Vetor Nulo
B0_R_OA = Vector(0)
# Vetor que liga os pontos A e B expresso no referencial móvel B2
B2_R_AB = L_1 * B2.x
# Vetor que liga os pontos B e C expresso no referencial óel B3
B3_R_BC = L_2 * B3.x

# Cinemática do ponto A em relação ao referencial B0
R_A = B0_R_OA
V_A = time_derivative(R_A, B0)
A_A = time_derivative(V_A, B0)

# Cinemática do ponto B em relação ao referencial B0
R_B = R_A + B2_R_AB.express(B0)
V_B = time_derivative(R_B, B0)
A_B = time_derivative(V_B, B0)

# Cinemática do ponto C em relação ao referencial B0
R_C = B3_R_BC.express(B0)
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
