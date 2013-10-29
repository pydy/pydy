#!/usr/bin/env python

# standard library
import time
import os
import shutil
import glob

# external libraries
from scipy.integrate import odeint
from numpy import hstack, ones, pi, linspace, array
import matplotlib.pyplot as plt
from pydy_code_gen.code import numeric_right_hand_side
from pydy_code_gen.tests.models import \
    generate_n_link_pendulum_on_cart_equations_of_motion

methods = ['lambdify', 'theano', 'cython']

derivation_times = {}
integration_times = {method: {} for method in methods}
code_generation_times = {method: {} for method in methods}

for n in range(1, 11):

    title = "Pendulum with {} links.".format(n)
    print(title)
    print('=' * len(title))

    start = time.time()
    results = \
        generate_n_link_pendulum_on_cart_equations_of_motion(n,
                                                             cart_force=False)
    derivation_times[n] = time.time() - start
    print('The derivation took {} seconds.\n'.format(derivation_times[n]))

    # Define the numerical values: parameters, time, and initial conditions
    arm_length = 1. / n
    bob_mass = 0.01 / n
    parameter_vals = [9.81, 0.01 / n]
    for i in range(n):
        parameter_vals += [arm_length, bob_mass]
    x0 = hstack((0, pi / 2 * ones(len(results[3]) - 1), 1e-3 *
                 ones(len(results[4]))))
    args = {'constants': array(parameter_vals),
            'num_coordinates': n + 1}
    t = linspace(0, 10, 1000)

    for method in methods:

        subtitle = "Generating with {} method.".format(method)
        print(subtitle)
        print('-' * len(subtitle))
        start = time.time()
        right_hand_side = numeric_right_hand_side(*results, generator=method)
        code_generation_times[method][n] = time.time() - start
        print('The code generation took {} seconds.'.format(
            code_generation_times[method][n]))

        start = time.time()
        y = odeint(right_hand_side, x0, t, args=(args,))
        integration_times[method][n] = time.time() - start
        print('ODE integration took {} seconds.\n'.format(
            integration_times[method][n]))

# clean up the cython crud
files = glob.glob('multibody_system*')
for f in files:
    os.remove(f)
shutil.rmtree('bulid')

# plot the results
fig, ax = plt.subplots(3, 1)

ax[0].plot(derivation_times.keys(), derivation_times.values(),
           label='Symbolic Derivation')
ax[0].set_title('Symbolic Derivation Time')

for k, v in code_generation_times.items():
    ax[1].plot(v.keys(), v.values(), label=k)
ax[1].set_title('Code Generation Time')
ax[1].legend()

for k, v in integration_times.items():
    ax[2].plot(v.keys(), v.values(), label=k)
ax[2].set_title('Integration Time')
ax[2].legend()

for a in ax.flatten():
    a.set_xlabel('Number of links.')
    a.set_ylabel('Time [s]')

fig.savefig('benchmark-results.png')
