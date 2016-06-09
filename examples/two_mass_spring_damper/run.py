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

mass_matrix = Matrix([[1, 0,  0,  0],
                      [0, 1,  0,  0],
                      [0, 0, m1,  0],
                      [0, 0,  0, m2]])
forcing = Matrix([u1,
                  u2,
                  -(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
                  k2*x1 - k2*x2 + c2*u1 - c2*u2])

# Form the remaining input lists for the equations of motion class

coordinates = (x1, x2)
speeds = (u1, u2)

# Initialize the equation of motion class

eom = eombase.EOM(coordinates, speeds, mass_matrix, forcing)

# Pass the equations of motion to system for simulations

sys = System(eom)

# Set up the system for simulation

sys.times = linspace(0, 10, num=100)
sys.constants = {m1: 10.0, c1: 10.0, k1: 10.0, m2: 5.0, c2: 5.0, k2: 5.0}
sys.initial_conditions = {x1: 1.0, x2: 0.0, u1: 0.0, u2: 0.0}

# Simulate the system

out = sys.integrate()

# Graph the results
plt.plot(sys.times, out[:, 1])  # Don't know output format of sys.integrate()
plt.plot(sys.times, out[:, 2])
plt.show()
