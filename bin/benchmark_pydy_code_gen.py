#!/usr/bin/env python

# standard library
import time
import os
import shutil
import glob

# external libraries
from scipy.integrate import odeint
from numpy import hstack, ones, pi, linspace, array, zeros, zeros_like
import matplotlib.pyplot as plt
from pydy.codegen.code import generate_ode_function
from pydy.codegen.tests.models import \
    generate_n_link_pendulum_on_cart_equations_of_motion


def run_benchmark(max_num_links, num_time_steps=1000):
    """Runs the n link pendulum derivation, code generation, and integration
    for each n up to the max number provided and generates a plot of the
    results."""

    methods = ['lambdify', 'theano', 'cython']

    link_numbers = range(1, max_num_links + 1)

    derivation_times = zeros(len(link_numbers))
    integration_times = zeros((max_num_links, len(methods)))
    code_generation_times = zeros_like(integration_times)

    for j, n in enumerate(link_numbers):

        title = "Pendulum with {} links.".format(n)
        print(title)
        print('=' * len(title))

        start = time.time()
        results = \
            generate_n_link_pendulum_on_cart_equations_of_motion(n, cart_force=False)
        derivation_times[j] = time.time() - start
        print('The derivation took {:1.5f} seconds.\n'.format(derivation_times[j]))

        # Define the numerical values: parameters, time, and initial conditions
        arm_length = 1. / n
        bob_mass = 0.01 / n
        parameter_vals = [9.81, 0.01 / n]
        for i in range(n):
            parameter_vals += [arm_length, bob_mass]

        # odeint arguments
        x0 = hstack((0, pi / 2 * ones(len(results[3]) - 1), 1e-3 *
                    ones(len(results[4]))))
        args = {'constants': array(parameter_vals)}
        t = linspace(0, 10, num_time_steps)

        for k, method in enumerate(methods):

            subtitle = "Generating with {} method.".format(method)
            print(subtitle)
            print('-' * len(subtitle))
            start = time.time()
            evaluate_ode = generate_ode_function(*results,
                                                      generator=method)
            code_generation_times[j, k] = time.time() - start
            print('The code generation took {:1.5f} seconds.'.format(
                code_generation_times[j, k]))

            start = time.time()
            odeint(evaluate_ode, x0, t, args=(args,))
            integration_times[j, k] = time.time() - start
            print('ODE integration took {:1.5f} seconds.\n'.format(
                integration_times[j, k]))

        del results, evaluate_ode

    # clean up the cython crud
    files = glob.glob('multibody_system*')
    for f in files:
        os.remove(f)
    shutil.rmtree('build')

    # plot the results
    fig, ax = plt.subplots(3, 1, sharex=True)

    ax[0].plot(link_numbers, derivation_times)
    ax[0].set_title('Symbolic Derivation Time')

    ax[1].plot(link_numbers, code_generation_times)
    ax[1].set_title('Code Generation Time')
    ax[1].legend(methods, loc=2)

    ax[2].plot(link_numbers, integration_times)
    ax[2].set_title('Integration Time')
    ax[2].legend(methods, loc=2)

    for a in ax.flatten():
        a.set_ylabel('Time [s]')

    ax[-1].set_xlabel('Number of links')

    plt.tight_layout()

    fig.savefig('benchmark-results.png')

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description='Run the n link pendulum benchmark.')

    parser.add_argument('max_num_links', type=int,
        help="The maximum number of links to compute.")

    parser.add_argument('num_time_steps', type=int,
        help="The number of integration time steps.")

    args = parser.parse_args()

    run_benchmark(args.max_num_links, args.num_time_steps)
