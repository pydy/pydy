#!/usr/bin/env python

"""This example shows how you can parallelize ODE integration of a generated
rhs function using the multiprocessing module. For example, this could be
useful in genetic algorithm optimization routines to parallelize the
evaluation of the cost function at each iteration."""

from multiprocessing import Pool

import numpy as np
from scipy.integrate import odeint
from pydy.models import multi_mass_spring_damper

print('Generating EoMs')
sys = multi_mass_spring_damper(10)

print('Defining numerical values')
x = np.random.random(len(sys.states))
t = np.linspace(0.0, 10.0, 100000)
pp = np.random.random((16, len(sys.constants_symbols)))

print('Generating rhs function')
rhs = sys.generate_ode_function(generator='cython')

print('Defining wrapper functions.')
# These wrappers are used to provide a single argmument to the function in
# the Pool.map call. There doesn't seem to be an easy way to pass in
# multiple arguments to the function that is being mapped.


def rhs2(p):
    return rhs(x, t[0], p)


def odeint2(p):
    return odeint(rhs, x, t, args=(p,))

p = Pool(4)

print('Running parallel rhs')
res1 = p.map(rhs2, [p_i for p_i in pp])

print('Running parallel odeint.')
res2 = p.map(odeint2, [p_i for p_i in pp])
