#!/usr/bin/env python

# standard library
import time

# external libraries
from numpy import hstack, ones, pi, linspace, array, zeros, zeros_like
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pydy.models import n_link_pendulum_on_cart
from sympy import symbols


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
        sys = n_link_pendulum_on_cart(n, cart_force=False)

        m = symbols('m:{}'.format(n + 1))
        l = symbols('l:{}'.format(n))
        g = symbols('g')

        derivation_times[j] = time.time() - start
        print('The derivation took {:1.5f} seconds.\n'.format(derivation_times[j]))

        # Define the numerical values: parameters, time, and initial conditions
        arm_length = 1. / n
        bob_mass = 0.01 / n
        parameter_vals = [9.81, 0.01 / n]
        for i in range(n):
            parameter_vals += [arm_length, bob_mass]

        times = linspace(0, 10, num_time_steps)
        sys.times = times

        x0 = hstack(
            (0,
             pi / 2 * ones(len(sys.coordinates) - 1),
             1e-3 * ones(len(sys.speeds))))
        sys.initial_conditions = dict(zip(sys.states, x0))

        constants = [g, m[0]]
        for i in range(n):
            constants += [l[i], m[i + 1]]

        sys.constants = dict(zip(constants, array(parameter_vals)))

        for k, method in enumerate(methods):

            subtitle = "Generating with {} method.".format(method)
            print(subtitle)
            print('-' * len(subtitle))
            start = time.time()
            sys.generate_ode_function(generator=method)
            code_generation_times[j, k] = time.time() - start
            print('The code generation took {:1.5f} seconds.'.format(
                code_generation_times[j, k]))

            start = time.time()
            sys.integrate()
            integration_times[j, k] = time.time() - start
            print('ODE integration took {:1.5f} seconds.\n'.format(
                integration_times[j, k]))

        del sys

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
