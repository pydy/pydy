
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
                          specifieds={self.specified_symbol: np.ones(1)},
                          constants=self.constant_map)

    def test_init(self):

        # Check defaults for most attributes.
        # -----------------------------------
        sys = System(self.kane)

        assert sys.rhs == None
        assert sys.eom_method is self.kane
        assert sys.ode_solver is odeint
        assert len(sys.specifieds) == 1
        assert sys.specifieds.keys()[0] is dynamicsymbols('F')
        testing.assert_allclose(sys.specifieds.values(), np.zeros(1))
        testing.assert_allclose(sys.initial_conditions, np.zeros(2))
        assert sys.constants == dict(zip(symbols('m, k, c, g'),
                                     [1.0, 1.0, 1.0, 1.0]))

        # Specify a bunch of attributes during construction.
        # --------------------------------------------------
        ic = {symbol('x'): 3.6, symbol('v'): 4.3}
        sys = System(self.kane,
                     ode_solver=odeint,
                     specifieds={self.specified_symbol: np.ones(1)},
                     initial_conditions=ic,
                     constants=self.constant_map)

        assert sys.method is self.kane
        assert sys.specified_symbols is dynamicsymbols('F')
        testing.assert_allclose(sys.specified_values, np.ones(1))
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(sys.initial_conditions.values(), ic.values())
        assert sys.constants == self.constant_map

    def test_coordinates(self):
        assert self.sys.coordinates == self.kane._q

    def test_speeds(self):
        assert self.sys.speeds = self.kane._u

    def test_states(self):
        assert self.sys.states = self.kane._q + self.kane._u

    def test_constants(self):

        # User-specified numerical values for constants.
        constants = {symbol('m'): 3.0, symbol('c'): 4.5}

        # The System should fill in the remaining constants with the defaults.
        constants_filled_with_remaining_defaults = {
                symbol('m'): 3.0, symbol('c'): 4.5,
                symbol('c'): 1.0, symbol('g'): 1.0}

        # Construct a new system with our constants at construction time.
        # ---------------------------------------------------------------
        sys = System(self.kane, constants=constants)

        # The System class adds in defaults for the ones we did not provide.
        assert sys.constants == constants_filled_with_remaining_defaults


        # Set constants after construction.
        # ---------------------------------
        sys = System(self.kane)

        # All constants have the default value.
        assert sys.constants == dict(zip(symbols('m, k, c, g'),
                                     [1.0, 1.0, 1.0, 1.0]))
        sys.constants = constants

        # The System fills in defaults.
        assert sys.constants == constants_filled_with_remaining_defaults


        # Using the property as a dict.
        # -----------------------------
        # TODO


        # Provide a constant that isn't actually a constant.
        # --------------------------------------------------
        testing.assert_raises(ValueError,
                sys.constants, {symbol('x'): 1.3})
        testing.assert_raises(ValueError,
                sys.constants, {symbol('F'): 1.8})

    def test_specifieds(self):

        # TODO use the n-link pendulum so that there are multiple specifieds.
            # TODO try various different groupings of the specifieds.
            # that is, {(a, b, c): np.ones(3), d: 5.6}


        # Using the property as a dict.
        # -----------------------------
        # TODO


        # The specified function returns the correct number of values.
        # ------------------------------------------------------------
        testing.assert_raises(ValueError,
                self.sys.specifieds,
                {self.specified_symbol: lambda x, t: np.ones(10)}
                )

        # The specified symbol must exist in the equations of motion and not
        # be a state.
        testing.assert_raises(ValueError,
                self.sys.specifieds, {symbol('x'): 5.1}
                )
        testing.assert_raises(ValueError,
                self.sys.specifieds, {symbol('m'): 5.1}
                )

    def test_ode_solver(self):

        assert self.sys.ode_solver == odeint
        self.sys.ode_solver = max
        assert self.sys.ode_solver is max

        # ode_solver must be a function
        # -----------------------------
        testing.assert_raises(ValueError,
                self.sys.ode_solver, 5)

    def test_initial_conditions(self):

        # Partially provided ic's.
        ic = {symbol('v'): 6.1}

        # The System should fill in the remaining ic's with the defaults.
        ic_filled_with_remaining_defaults = {
                symbol('v'): 6.1, symbol('x'): 0.0}


        # Using the constructor.
        # ----------------------
        sys = System(self.kane, initial_conditions=ic)
        assert sys.initial_conditions == ic_filled_with_remaining_defaults


        # Set the attribute.
        # ------------------
        self.sys = ic
        assert self.sys.initial_conditions == ic_filled_with_remaining_defaults

        # The ordering must be correct. The dict might mess this up.
        # ----------------------------------------------------------
        # TODO


        # Using the property as a dict.
        # -----------------------------
        # TODO


        # Keys must be coords or speeds.
        # ------------------------------
        testing.assert_raises(ValueError,
                self.sys.initial_conditions={symbol('k'): 0.4}
                )
        testing.assert_raises(ValueError,
                self.sys.initial_conditions={symbol('F'): 7.3}
                )

    def test_generate_ode_function(self):

        # TODO
        rhs = self.sys.generate_ode_function()

        assert rhs is self.sys.evaluate_ode_function

        args = {'constants': self.constants.values(),
                'specified': np.zeros(1)}

        actual = rhs(np.ones(2), 0.0, args=(args,))

        testing.assert_allclose(actual,
                                np.array([0, 0]))
        raise Exception()

    def test_integrate(self):
        # TODO
        # TODO ensure that the reordering of initial conditions happens
        # properly.
        raise Exception()
