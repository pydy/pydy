
import numpy as np
from numpy import testing
from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols
from scipy.integrate import odeint

from ..system import System
from ..codegen.tests.models import \
    generate_mass_spring_damper_equations_of_motion as mass_spring_damper


class TestSystem():

    def setup(self):

        self.kane = mass_spring_damper(kane=True)
        self.specified_symbol = dynamicsymbols('F')
        self.constant_map = dict(zip(symbols('m, k, c, g'),
                                     [2.0, 1.5, 0.5, 9.8]))
        self.sys = System(self.kane,
                          code_gen_backend='lambdify',
                          specifieds={self.specified_symbol: np.ones(1)},
                          initial_conditions=np.zeros(2),
                          constants=self.constant_map)

    def test_init(self):

        sys = System(self.kane)

        assert sys.rhs == None
        assert sys.method is self.kane
        assert sys.code_gen_backend == 'lambdify'
        assert sys.ode_solver is odeint
        assert len(sys.specifieds) == 1
        assert sys.specifieds.keys()[0] is dynamicsymbols('F')
        testing.assert_allclose(sys.specifieds.values(), np.zeros(1))
        testing.assert_allclose(sys.initial_conditions, np.zeros(2))
        assert sys.constants == dict(zip(symbols('m, k, c, g'),
                                     [1.0, 1.0, 1.0, 1.0]))

        sys = System(self.kane,
                     code_gen_backend='lambdify',
                     ode_solver=odeint,
                     specified_symbols=(self.specified_symbol,),
                     specified_values=np.ones(1),
                     initial_conditions=np.zeros(2),
                     constants=self.constant_map)

        assert sys.method is self.kane
        assert sys.code_gen_backend == 'lambdify'
        assert sys.specified_symbols is dynamicsymbols('F')
        testing.assert_allclose(sys.specified_values, np.ones(1))
        testing.assert_allclose(sys.initial_conditions, np.zeros(2))
        assert sys.constants == self.constant_map

    def test_constants(self):

        constants_map = {symbol('m'): 3.0, symbol('c'): 4.5}
        sys = System(self.kane, constants=constants_map)
        assert sys.constants_map == constants_map


        sys.backend == 'cython'

        # You must set both specified symbols and values.
        assert_raises(ValueError,
                      System(self.kane,
                             specified_symbols=(self.specified_symbol)))

        assert_raises(ValueError,
                      System(self.kane,
                             specified_values=(self.specified_symbol)))

        # Check to make sure the specified function returns the correct
        # number of values.
        assert_raises(ValueError,
                      System(self.kane,
                             specified_symbols=(self.specified_symbol),
                             specified_values=lambda x, t: np.ones(10))


        # The specified symbol must exist in the equations of motion and not
        # be a state.
        assert_raises(ValueError,
                      System(self.kane,
                             specified_symbols=(dynamicsymbols('G')),
                             specified_values=np.ones(1)))

        # The same number of specified symbols and values are needed.
        assert_raises(ValueError,
                      System(self.kane,
                             specified_symbols=(self.specified_symbol),
                             specified_values=np.ones(2))


    def test_generate_ode_function(self):

        rhs = self.sys.generate_ode_function()

        assert rhs is self.rhs

        args = {'constants': self.constant_map.values,
                'specified': np.zeros(1)}

        actual = rhs(np.ones(2), 0.0, args=(args,))

        testing.assert_allclose(actual,
                                np.array([ , ]))

    def find_constants(self):

        constant_dict = self.sys._find_constants()

        assert constant_dict = self.constant_map

    def find_specified_symbols(self):

        specified_dict = self.sys._find_specified_symbols()

        assert specified_dict = {self.dynamic_symbol: 0.0}

    def test_integrate(self):

    """
    method : sympy.physics.mechanics.KanesMethod or
                sympy.physics.mechanics.LagrangesMethod
        The method used to generate the equations of motion.
    constants : dict, optional (default: all 1.0)
    specifieds : dict, optional (default: all 0)
    code_gen_backend : str, optional (default: 'lambdify')
    ode_solver : function, optional (default: scipy.integrate.odeint)
    initial_conditions : dict, optional (default: all zero)
        """
