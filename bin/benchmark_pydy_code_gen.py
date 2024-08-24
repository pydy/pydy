#!/usr/bin/env python

# standard library
import time

# external libraries
from numpy import hstack, ones, pi, linspace, array, zeros, zeros_like, nan
from pydy.models import n_link_pendulum_on_cart
from sympy import symbols
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')


def run_benchmark(max_num_links, num_time_steps=1000, duration=10.0):
    """Runs the n link pendulum derivation, code generation, and integration
    for each n up to the max number provided and generates a plot of the
    results."""

    methods = ['lambdify', 'cython', 'theano']

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
        msg = 'The derivation took {:1.5f} seconds.\n'
        print(msg.format(derivation_times[j]))

        # Define the numerical values: parameters, time, and initial conditions
        arm_length = 1.0 / n
        bob_mass = 0.01 / n
        parameter_vals = [9.81, 0.01 / n]
        for i in range(n):
            parameter_vals += [arm_length, bob_mass]

        times = linspace(0.0, duration, num=num_time_steps)
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
            try:
                sys.generate_ode_function(generator=method, cse=True)
            # ImportError: Theano or Cython not installed
            # AttributeError: Theano doesn't work with new NumPy versions
            except (ImportError, AttributeError):
                print("Skipped {} due to error.\n".format(method))
                code_generation_times[j, k] = nan
                integration_times[j, k] = nan
            else:
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
    fig, ax = plt.subplots(3, 1, sharex=True, layout='constrained')

    ax[0].plot(link_numbers, derivation_times)
    ax[0].set_title('Symbolic Derivation Time')

    ax[1].plot(link_numbers, code_generation_times)
    ax[1].set_yscale('log')
    ax[1].set_title('Code Generation Time')
    ax[1].legend(methods, loc=2)

    ax[2].plot(link_numbers, integration_times)
    ax[2].set_yscale('log')
    ax[2].set_title('Integration Time')
    ax[2].legend(methods, loc=2)

    for a in ax.flatten():
        a.set_ylabel('Time [s]')

    ax[-1].set_xlabel('Number of links')

    fig.savefig('benchmark-results.png')


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description='Run the n link pendulum benchmark.')

    parser.add_argument('max_num_links', type=int,
                        help="The maximum number of links to compute.")

    parser.add_argument('num_time_steps', type=int,
                        help="The number of integration time steps.")

    parser.add_argument('duration', type=float,
                        help="The simulation duration.")

    args = parser.parse_args()

    run_benchmark(args.max_num_links, args.num_time_steps, args.duration)
