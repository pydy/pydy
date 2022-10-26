#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
from pydy.codegen.code import generate_ode_function

from derive import *

# List the symbolic arguments
# ===========================

# Constants
# ---------

constants = {r: 5.0, m: 10.0, g: 9.81}

# Time-varying
# ------------

coordinates = [q1, q2, q3, q4, q5]
speeds = [u1, u2, u3]
states = coordinates + speeds

xdot_function = generate_ode_function(M, F,
                                      constants.keys(),
                                      coordinates, speeds)

initial_conditions = {q1: 0.0,
                      q2: np.deg2rad(5.0),
                      q3: 0.0,
                      q4: 0.0,
                      q5: 0.0,
                      u1: 0.0,
                      u2: 10.0 / constants[r],
                      u3: 0.0}

x0 = [initial_conditions[s] for s in states]

args = {'constants': constants}

# Simulate
# ========

frames_per_sec = 60
final_time = 10.0

t = np.linspace(0.0, final_time, num=int(final_time * frames_per_sec) + 1)
x = odeint(xdot_function, x0, t, args=(args,))
