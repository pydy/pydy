#!/usr/bin/env python

"""This example shows how you can parallelize ODE integration of a generated
ODE function using the multiprocessing module. The example shows how you can
both evaluate the right hand side function and integrate the equations of
motion with different model parameters while spreading the independent
computations over the number of CPUs on the computer. For example, this
could be useful in genetic algorithm optimization routines to parallelize
the evaluation of the cost function at each iteration."""

from multiprocessing import Pool

import numpy as np
from scipy.integrate import odeint
from pydy.models import multi_mass_spring_damper

print('Generating equations of motion')
sys = multi_mass_spring_damper(10)

print('Defining numerical values')
x = np.random.random(len(sys.states))
t = np.linspace(0.0, 10.0, 100000)
# 16 different parameter sets to evaluate in parallel.
p_set = np.random.random((16, len(sys.constants_symbols)))

print('Generating the ODE function')
rhs = sys.generate_ode_function(generator='cython')

print('Defining wrapper functions')
# These wrappers are used to provide a single argmument to the function in
# the Pool.map() call. There doesn't seem to be an easy way to pass in
# multiple arguments to the function that is being mapped.


def rhs_wrapper(p):
    return rhs(x, t[0], p)


def odeint_wrapper(p):
    return odeint(rhs, x, t, args=(p,))

pool = Pool()

print('Running rhs evalutions in parallel')
res1 = pool.map(rhs_wrapper, [p for p in p_set])

print('Running odeint evaluations in parallel')
res2 = pool.map(odeint_wrapper, [p for p in p_set])
