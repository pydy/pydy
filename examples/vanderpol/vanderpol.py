import numpy as np, sympy as sp, matplotlib.pyplot as plt, numpy.matlib
from pydy.system import System
from sympy.physics.mechanics import dynamicsymbols

x = sp.Matrix(dynamicsymbols('x1:3'))
x1, x2 = x

e = sp.symbols('epsilon')

sys = System(x, sp.Matrix([x2, -x1+e*(1-x1**2)*x2]))

sys.constants = {e: 5}
sys.initial_conditions = {x1: 1, x2: 1}
sys.times = np.linspace(0, 10, 500)
res = sys.integrate()

plt.plot(sys.times,res)
plt.legend([sp.latex(s, mode='inline') for s in sys.states])