#!/usr/bin/env python

"""This example shows how you can parallelize ODE integration of a generated
rhs function using joblib. For example, this could be useful in genetic
algorithm optimization routines to parallelize the evaluation of the cost
function at each iteration."""

import numpy as np
from scipy.integrate import odeint
from joblib import Parallel, delayed
from pydy.models import n_link_pendulum_on_cart

print('Generating EoMs')
sys = n_link_pendulum_on_cart(10, False, False)

print('Defining numerical values')
x = np.random.random(len(sys.states))
t = np.linspace(0.0, 10.0, 100000)
pp = np.random.random((16, len(sys.constants_symbols)))

print('Generating rhs function')
rhs = sys.generate_ode_function(generator='cython')

print('Defining wrapper functions.')
# These two wrapper functions seem to be required for things to pickle in
# Joblib and I'm not sure why. This is not the case if the rhs function was
# defined as a normal Python function instead of being generated with
# lambdify or cython.


def rhs2(p):
    return rhs(x, t[0], p)


def odeint2(p):
    return odeint(rhs, x, t, args=(p,))


print('Running parallel rhs')
res1 = Parallel(n_jobs=-1)(delayed(rhs2)(p_i) for p_i in pp)

print('Running parallel odeint.')
res2 = Parallel(n_jobs=-1)(delayed(odeint2)(p_i) for p_i in pp)
