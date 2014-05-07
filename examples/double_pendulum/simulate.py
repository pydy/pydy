"""
This file will use pydy.codegen to simulate the double pendulum.

"""

from numpy import concatenate, array, linspace
from pydy.codegen.code import generate_ode_function
from scipy.integrate import odeint

from double_pendulum import *

# List the symbolic arguments
# ===========================

# Constants
# ---------

constants = {l: 10.0, m: 10.0, g: 9.81}

# Time-varying
# ------------

coordinates = [q1, q2]

speeds = [u1, u2]


# Generate function that returns state derivatives
# ================================================

xdot_function = generate_ode_function(mass_matrix, forcing_vector,
        constants.keys(), coordinates, speeds)


# Specify numerical quantities
# ============================

initial_coordinates = [1.0, 0.0]
initial_speeds = [0.0, 0.0]
x0 = concatenate((initial_coordinates, initial_speeds), axis=1)

args = {'constants': constants.values()}


# Simulate
# ========

frames_per_sec = 60
final_time = 5.0

t = linspace(0.0, final_time, final_time * frames_per_sec)
x = odeint(xdot_function, x0, t, args=(args,))
