import time
from scipy.integrate import odeint
from numpy import hstack, ones, pi, linspace

from benchmark_functions import generate_equations, numeric_right_hand_side

for n in range(1, 11):

    print("n = {}".format(n))

    start = time.time()
    kane, parameters = generate_equations(n)
    print('The derivation took {} seconds.'.format(time.time() - start))

    # Define the numerical values: paramters, time, and initial conditions
    arm_length = 1. / n
    bob_mass = 0.01 / n
    parameter_vals = [9.81, 0.01 / n]
    for i in range(n):
        parameter_vals += [arm_length, bob_mass]
    x0 = hstack(( 0, pi / 2 * ones(len(kane._q) - 1), 1e-3 * ones(len(kane._u)) ))
    t = linspace(0, 10, 1000)

    methods = ['lambdify', 'theano', 'autowrap']
    for method in methods:
        print("Running with {} method.".format(method))
        # TODO : time the code generation part
        # TODO : it'd also be interesting to know the time it takes to call solve
        # as solve and the integration can be pushed to the faster languages too.
        # We need to have a way to seperate that if the odes are pass in before the
        # matrix is solved.
        right_hand_side = numeric_right_hand_side(kane, parameters, generator=method)
        #print('Starting odeint')
        #start = time.time()
        #y = odeint(right_hand_side, x0, t, args=(parameter_vals,))
        #print('ODE integration took {} seconds.'.format(time.time() - start))
