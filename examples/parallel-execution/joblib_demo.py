#!/usr/bin/env python

"""This example shows how you can parallelize ODE integration of a generated
ODE function using joblib. The example shows how you can both evaluate the
right hand side function and integrate the equations of motion with
different model parameters while spreading the independent computations over
the number of CPUs on the computer. For example, this could be useful in
genetic algorithm optimization routines to parallelize the evaluation of the
cost function at each iteration."""

import numpy as np
from scipy.integrate import odeint
from joblib import Parallel, delayed
from pydy.models import n_link_pendulum_on_cart

print('Generating equations of motion')
sys = n_link_pendulum_on_cart(10, False, False)

print('Defining numerical values')
x = np.random.random(len(sys.states))
t = np.linspace(0.0, 10.0, 100000)
p_set = np.random.random((16, len(sys.constants_symbols)))

print('Generating the ODE function')
rhs = sys.generate_ode_function(generator='cython')

print('Defining wrapper functions')
# These two wrapper functions seem to be required for things to pickle in
# Joblib and I'm not sure why. This is not the case if the rhs function was
# defined as a normal Python function instead of being generated with
# lambdify or cython.


def rhs_wrapper(p):
    return rhs(x, t[0], p)


def odeint_wrapper(p):
    return odeint(rhs, x, t, args=(p,))

# The following calls to Parallel and delayed are recommended by the joblib
# documentation. See the documentation for more information.

print('Running rhs evalutions in parallel')
res1 = Parallel(n_jobs=-1)(delayed(rhs_wrapper)(p) for p in p_set)

print('Running odeint evaluations in parallel')
res2 = Parallel(n_jobs=-1)(delayed(odeint_wrapper)(p) for p in p_set)
