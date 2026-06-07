# This code will go over the manual input of the equations of motion for the
# simple pendulum into eombase using x and y coordinates instead of theta

# The equations of motion are formed at
# http://nbviewer.jupyter.org/github/bmcage/odes/blob/master/docs/ipython/Planar%20Pendulum%20as%20DAE.ipynb
from sympy import symbols, Matrix
from sympy.physics.mechanics import dynamicsymbols, eombase

# Define the dynamic symbols
x, y, u, v, lam = dynamicsymbols('x y u v lambda')

# Define the constant symbols
m, l, g = symbols('m l g')

# Define the mass matrix and forcing vector
mm = Matrix([[1, 0, -x/m],
             [0, 1, -y/m],
             [0, 0, l**2/m]])
f = Matrix([0, 0, u**2 + v**2 - g*y])

mm_full = Matrix([[1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0],
                  [0, 0, 1, 0, -x/m],
                  [0, 0, 0, 1, -y/m],
                  [0, 0, 0, 0, l**2/m]])
f_full = Matrix([u, v, 0, 0, u**2 + v**2 - g*y])

# Form the rhs of q' = G(q, u, t, r, p). Kinematic equation
G = Matrix([u, v])

# Form the rhs of the dynamics equations
RHS = mm_full.inv()*f_full

# Create a list specifing the rows conatining algebraic rather than differential
# constraints
alg_con = [2]
alg_con_full = [4]

# Create the interable of states, coordinates and speeds
states = (x, y, u, v, lam)

# Initialize the equation of motion class using the three forms accepted by
# ODEFunctionGenerator
#    [1] x' = F(x, t, r, p)
#
#    [2] M(x, p) x' = F(x, t, r, p)
#
#    [3] M(q, p) u' = F(q, u, t, r, p)
#        q' = G(q, u, t, r, p)

eom1 = eombase.EOM(states=states, rhs=RHS, alg_con=alg_con_full)
eom2 = eombase.EOM(states=states, mass_matrix_full=mm_full, forcing_full=f_full,
                   alg_con=alg_con_full)
eom3 = eombase.EOM(states=states, mass_matrix=mm, forcing=f, kinematics=G,
                   alg_con=alg_con)
