"""
This is a script to test the time taken by the generated
rhs function in ode_function_generators.

We allow the user to specify what type, the extra args
should be for both the constants and the specifieds. The constants
can be None, 'array', or 'dictionary'. The specifieds can be None,
'array', 'function', or 'dictionary'. Thus, we have 12 permutations.

The script outputs the average time taken to numerically
evaluate the right hand side of the first order differential
equation for each of the 12 permutations.

"""

import timeit

import numpy as np
import sympy as sm

from pydy import models
from pydy.codegen.ode_function_generators import LambdifyODEFunctionGenerator

sys = models.n_link_pendulum_on_cart(3, True, True)

right_hand_side = sys.eom_method.rhs()

constants = list(sm.ordered(sys.constants_symbols))

specifieds = list(sm.ordered(sys.specifieds_symbols))

constants_arg_types = ['array', 'dictionary']
specifieds_arg_types = ['array', 'function', 'dictionary']

p_array = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
p_dct = dict(zip(constants, p_array))

p = {}

p['array'] = p_array
p['dictionary'] = p_dct

r_array = np.array([1.0, 2.0, 3.0, 4.0])
# r_dct_1 = dict(zip(specifieds, r_array))
# r_dct_2 = {tuple(specifieds):
#           lambda x, t: r_array}
r_dct_3 = {specifieds[0]: lambda x, t: np.ones(1),
           (specifieds[3], specifieds[1]):
           lambda x, t: np.array([4.0, 2.0]),
           specifieds[2]: 3.0 * np.ones(1)}
r_func = lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])

r = {}

r['array'] = r_array
r['dictionary'] = r_dct_3
r['function'] = r_func

x = np.array([0.125, 0.250, 0.375, 0.456,
              0.123, 0.457, 0.999, 0.192])
itr = 1000
print('The time taken by rhs function in {} iterations'.format(itr))

for p_arg_type in constants_arg_types:
    for r_arg_type in specifieds_arg_types:

        g = LambdifyODEFunctionGenerator(right_hand_side,
                                         sys.coordinates,
                                         sys.speeds,
                                         constants,
                                         specifieds=specifieds,
                                         constants_arg_type=p_arg_type,
                                         specifieds_arg_type=r_arg_type)

        rhs = g.generate()

        time = timeit.repeat("rhs(x, 0.0, r[r_arg_type], p[p_arg_type])",
                             "from __main__ import rhs,x,r,p,"
                             "r_arg_type,p_arg_type",
                             number=itr)

        print('For constants argument type - "{p_arg}" and '
              'specifieds argument type - "{r_arg}" is '
              .format(p_arg=p_arg_type, r_arg=r_arg_type))

        print(sum(time)/3)
