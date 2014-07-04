
from numpy import concatenate, array, linspace, pi
from scipy.integrate import odeint
from pydy.codegen.code import generate_ode_function

from derive import *

# List the symbolic arguments
# ===========================

# Constants
# ---------

constants = {r: 10*0.25 * 3.141459, m: 10.0, g: 9.81}

# Time-varying
# ------------

coordinates = [q1, q2, q3, q4, q5]
speeds = [u1, u2, u3]

xdot_function = generate_ode_function(M, F, constants.keys(), coordinates, speeds)

initial_coordinates = [0.0, 0.01, 0.0, 0.0, 0.0]
initial_speeds = [0.0, 5.0*0.25*pi, 0.0]
x0 = concatenate((initial_coordinates, initial_speeds), axis=1)

args = {'constants': constants.values()}

# Simulate
# ========

frames_per_sec = 60
final_time = 10.0

t = linspace(0.0, final_time, final_time * frames_per_sec)
x = odeint(xdot_function, x0, t, args=(args,))
