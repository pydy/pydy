import matplotlib.pyplot as plt
from numpy import linspace
from pydy.system import System
from sympy import symbols, Matrix
from sympy.physics.mechanics import dynamicsymbols, eombase

# Define the dynamic symbols

x1, x2, u1, u2 = dynamicsymbols('x1 x2 u1 u2')

# Define the constant symbols

m1, c1, k1 = symbols('m1 c1 k1')
m2, c2, k2 = symbols('m2 c2 k2')

# Define the mass matrix and forcing vector

mm = Matrix([[m1,  0],
             [0,  m2]])
f = Matrix([-(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
            k2*x1 - k2*x2 + c2*u1 - c2*u2])

mm_full = Matrix([[1, 0,  0,  0],
                  [0, 1,  0,  0],
                  [0, 0, m1,  0],
                  [0, 0,  0, m2]])
f_full = Matrix([u1,
                 u2,
                 -(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
                 k2*x1 - k2*x2 + c2*u1 - c2*u2])

# Form the rhs of q' = G(q, u, t, r, p). Kinematic equation

G = Matrix([u1, u2])

# Form the rhs of the dynamics equations

RHS = mm_full.inv()*f_full

# Form the remaining input lists for the equations of motion class

coordinates = (x1, x2)
speeds = (u1, u2)
states = (x1, x2, u1, u2)

# Form the kinetic energy to attach to one of the EOM instances as an output
# equation

KE = 0.5 * (m1*u1**2 + m2*u2**2)
out_eqns = {"kinetic_energy": KE}

# Initialize the equation of motion class using the three forms accepted by
# ODEFunctionGenerator
#    [1] x' = F(x, t, r, p)
#
#    [2] M(x, p) x' = F(x, t, r, p)
#
#    [3] M(q, p) u' = F(q, u, t, r, p)
#        q' = G(q, u, t, r, p)

eom1 = eombase.EOM(coordinates=coordinates, speeds=speeds, rhs=RHS,
                   output_eqns=out_eqns)
eom2 = eombase.EOM(states=states, mass_matrix_full=mm_full, forcing_full=f_full,
                   num_coordinates=2)
eom3 = eombase.EOM(states=states, mass_matrix=mm, forcing=f, kinematics=G)

# Pass the equations of motion to system for simulations

sys = System(eom2)

# Set up the system for simulation

sys.times = linspace(0, 10, num=100)
sys.constants = {m1: 10.0, c1: 10.0, k1: 10.0, m2: 5.0, c2: 5.0, k2: 5.0}
sys.initial_conditions = {x1: 1.0, x2: 0.0, u1: 0.0, u2: 0.0}

# Simulate the system

out = sys.integrate()


# Display the kinetic energy from the first EOM class
# Returns a dictionary whith the same keys as output_eqns but the values are the
# numerical results of the simulation

print(eom1.output_eqns_results)

# Graph the results
plt.plot(sys.times, out[:, 1])  # Don't know output format of sys.integrate()
plt.plot(sys.times, out[:, 2])
plt.show()

# Showcase System's ability to plot the coordinates
sys.plot_coordinates()
