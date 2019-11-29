#!/usr/bin/env python

# external
from numpy import radians, linspace, hstack, zeros, ones
from scipy.integrate import odeint
from pydy.codegen.ode_function_generators import generate_ode_function

# local
from derive import l, m_bob, m_link, Ixx, Iyy, Izz, g, kane, q, u

param_syms = []
for par_seq in [l, m_bob, m_link, Ixx, Iyy, Izz, (g,)]:
    param_syms += list(par_seq)

# All of the links and bobs will have the same numerical values for the
# parameters.

link_length = 10.0  # meters
link_mass = 10.0  # kg
link_radius = 0.5  # meters
link_ixx = 1.0 / 12.0 * link_mass * (3.0 * link_radius**2 + link_length**2)
link_iyy = link_mass * link_radius**2
link_izz = link_ixx

particle_mass = 5.0  # kg
particle_radius = 1.0  # meters

# Create a list of the numerical values which have the same order as the
# list of symbolic parameters.
param_vals = [link_length for x in l] + \
             [particle_mass for x in m_bob] + \
             [link_mass for x in m_link] + \
             [link_ixx for x in list(Ixx)] + \
             [link_iyy for x in list(Iyy)] + \
             [link_izz for x in list(Izz)] + \
             [9.8]

# A function that evaluates the right hand side of the set of first order
# ODEs can be generated.
print("Generating numeric right hand side.")
right_hand_side = generate_ode_function(kane.forcing_full, q, u, param_syms,
                                        mass_matrix=kane.mass_matrix_full,
                                        generator='cython')

# To simulate the system, a time vector and initial conditions for the
# system's states is required.
t = linspace(0.0, 60.0, num=600)
x0 = hstack((ones(6) * radians(10.0), zeros(6)))

print("Integrating equations of motion.")
state_trajectories = odeint(right_hand_side, x0, t, args=(dict(zip(param_syms,
                                                                   param_vals)),
                                                          ))
print("Integration done.")
