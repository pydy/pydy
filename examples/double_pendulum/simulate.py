"""
This file will use pydy.codegen to simulate the double pendulum.

"""

from numpy import linspace
from pydy.system import System

from double_pendulum import *

constants = {l: 10.0, m: 10.0, g: 9.81}

initial_conditions = {q1: 1.0, q2: 0.0, u1: 0.0, u2: 0.0}

sys = System(KM, constants=constants,
        initial_conditions=initial_conditions)

frames_per_sec = 60
final_time = 5.0

times = linspace(0.0, final_time, final_time * frames_per_sec)
sys.times = times
x = sys.integrate()
