#!/usr/bin/env python

# standard library
import time
import timeit

# external libraries
from pydy.models import n_link_pendulum_on_cart
from sympy import symbols
import matplotlib.pyplot as plt
import numpy as np


def run_benchmark(max_num_links, num_time_steps=1000, duration=10.0):
    """Runs the n link pendulum derivation, code generation, and integration
    for each n up to the max number provided and generates a plot of the
    results."""

    methods = ['lambdify', 'cython', 'theano', 'symjit', 'cython:sympy']

    link_numbers = range(1, max_num_links + 1)

    derivation_times = np.zeros(len(link_numbers))
    rhs_times = np.zeros((max_num_links, len(methods)))
    integration_times = np.zeros_like(rhs_times)
    code_generation_times = np.zeros_like(rhs_times)

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
        arm_length = 1.0/n
        bob_mass = 0.01/n
        parameter_vals = [9.81, 0.01/n]
        for i in range(n):
            parameter_vals += [arm_length, bob_mass]

        times = np.linspace(0.0, duration, num=num_time_steps)
        sys.times = times

        x0 = np.hstack((
            0.0,
            np.pi/2*np.ones(len(sys.coordinates) - 1),
            1e-3*np.ones(len(sys.speeds)),
        ))
        sys.initial_conditions = dict(zip(sys.states, x0))

        constants = [g, m[0]]
        for i in range(n):
            constants += [l[i], m[i + 1]]

        sys.constants = dict(zip(constants, np.array(parameter_vals)))

        for k, method in enumerate(methods):

            subtitle = "Generating with {} method.".format(method)
            print(subtitle)
            print('-' * len(subtitle))
            start = time.time()
            try:
                if method == 'cython:sympy':
                    rhs = sys.generate_ode_function(generator='cython',
                                                    linear_sys_solver='sympy',
                                                    cse=True)
                else:
                    rhs = sys.generate_ode_function(generator=method, cse=True)
            # ImportError: Theano or Cython not installed
            # AttributeError: Theano doesn't work with new NumPy versions
            except (ImportError, AttributeError):
                print("Skipped {} due to error.\n".format(method))
                code_generation_times[j, k] = np.nan
                integration_times[j, k] = np.nan
            else:
                code_generation_times[j, k] = time.time() - start
                print('The code generation took {:1.5f} seconds.'.format(
                    code_generation_times[j, k]))

                p_vals = np.array(parameter_vals)
                rhs_time = timeit.timeit(lambda: rhs(x0, 0.1, p_vals),
                                         number=1000)
                rhs_times[j, k] = rhs_time
                print('rhs() evaluation took {:1.5f} seconds.'.format(
                    rhs_time))

                start = time.time()
                sys.integrate()
                integration_times[j, k] = time.time() - start
                print('ODE integration took {:1.5f} seconds.\n'.format(
                    integration_times[j, k]))

        del sys

    # plot the results
    fig, ax = plt.subplots(4, 1, sharex=True, layout='constrained',
                           figsize=(6, 6))

    ax[0].plot(link_numbers, derivation_times)
    ax[0].set_title('Symbolic Derivation Time')

    ax[1].plot(link_numbers, code_generation_times)
    ax[1].set_yscale('log')
    ax[1].set_title('Code Generation Time')
    ax[1].legend(methods, loc=2)

    ax[2].plot(link_numbers, rhs_times)
    ax[2].set_yscale('log')
    ax[2].set_title('ODE Evaluation Time')
    ax[2].legend(methods, loc=2)

    ax[3].plot(link_numbers, integration_times)
    ax[3].set_yscale('log')
    ax[3].set_title('Integration Time')
    ax[3].legend(methods, loc=2)

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
