
import numpy as np
from numpy import testing
from sympy import symbols
from sympy.physics.mechanics import dynamicsymbols
from scipy.integrate import odeint

from ..system import System
from ..codegen.tests.models import \
    generate_mass_spring_damper_equations_of_motion as mass_spring_damper
from ..codegen.tests.models import \
    generate_n_link_pendulum_on_cart_equations_of_motion as n_link_cart



class TestSystem():

    def setup(self):

        self.kane = mass_spring_damper(only_return_kane=True)
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
        ic = {dynamicsymbols('x'): 3.6, dynamicsymbols('v'): 4.3}
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
        assert self.sys.speeds == self.kane._u

    def test_states(self):
        assert self.sys.states == self.kane._q + self.kane._u

    def test_constants(self):

        # User-specified numerical values for constants.
        constants = {symbols('m'): 3.0, symbols('c'): 4.5}

        # The System should fill in the remaining constants with the defaults.
        constants_filled_with_remaining_defaults = {
                symbols('m'): 3.0, symbols('c'): 4.5,
                symbols('k'): 1.0, symbols('g'): 1.0}

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
        # Modifying hte dict directly does change the dict.
        sys.constants[symbols('m')] = 5.1
        assert sys.constants == dict(zip(symbols('m, k, c, g'),
                                     [5.1, 1.0, 1.0, 1.0]))
        # Putting in a non-constant key does not raise exception.
        sys.constants[dynamicsymbols('v')] = 3.5
        # Then, if we integrate, we do error-checking and we get an exception.
        testing.assert_raises(ValueError, sys.integrate, [0.0, 1.0])


        # Provide a constant that isn't actually a constant.
        # --------------------------------------------------
        testing.assert_raises(ValueError,
                sys.constants, {dynamicsymbols('x'): 1.3})
        testing.assert_raises(ValueError,
                sys.constants, {dynamicsymbols('F'): 1.8})

    def test_specifieds(self):

        sys = System(self.kane)
        assert sys.specifieds == {dynamicsymbols('F'), 0.0}

        # Using the property as a dict.
        # -----------------------------
        # Modifying hte dict directly does change the dict.
        sys.constants[dynamicsymbols('F')] = 5.1
        assert sys.constants == {dynamicsymbols('F'): 5.1}
        # Putting in a non-specified key does not raise exception.
        sys.constants[dynamicsymbols('v')] = 3.5
        # Then, if we integrate, we do error-checking and we get an exception.
        testing.assert_raises(ValueError, sys.integrate, [0.0, 1.0])

        # Putting in a value of the wrong length does not raise exception.
        sys.constants[dynamicsymbols('F')] = 3.1 * np.ones(2)
        # Then, if we integrate, we do error-checking and we get an exception.
        testing.assert_raises(ValueError, sys.integrate, [0.0, 1.0])


        # The specified function returns the correct number of values.
        # ------------------------------------------------------------
        testing.assert_raises(ValueError,
                self.sys.specifieds,
                {self.specified_symbol: lambda x, t: np.ones(10)}
                )

        # The specified symbol must exist in the equations of motion and not
        # be a state.
        testing.assert_raises(ValueError,
                self.sys.specifieds, {dynamicsymbols('x'): 5.1}
                )
        testing.assert_raises(ValueError,
                self.sys.specifieds, {symbols('m'): 5.1}
                )

        # Complex error-checking when using property as a dict.
        # -----------------------------------------------------
        km = n_link_cart(3, cart_force=True, joint_torques=True,
                only_return_kane=True)
        sys = System(km)
        spec_syms = sys.specifieds.keys()
        times = np.linspace(0, 0.5, 10)
        sys.specifieds = {spec_syms[0]: lambda x, t: np.ones(t),
                (spec_syms[3], spec_syms[1]): lambda x, t: np.array([4, 2]),
                spec_syms[2]: 3.0 * np.ones(1)}
        # These won't throw an exception.
        sys.specifieds[spec_syms[1]] = 7.1
        testing.assert_raises(ValueError, sys.integrate, times)

        sys.specifieds[(spec_syms[0], spec_syms[3])] = 5.8

        km = n_link_cart(3, cart_force=True, joint_torques=True,
                only_return_kane=True)
        sys = System(km)
        # This puts too many entries in the dict.
        sys.specifieds[(spec_syms[0], spec_syms[3])] = 5.8
        testing.assert_raises(ValueError, sys.integrate, times)

        # This gets rid of the previous default entries, and should work
        # properly.
        sys.specifieds.pop(spec_symbols[0])
        sys.specifieds.pop(spec_symbols[3])
        sys.integrate(times)


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
        ic = {dynamicsymbols('v'): 6.1}

        # The System should fill in the remaining ic's with the defaults.
        ic_filled_with_remaining_defaults = {
                dynamicsymbols('v'): 6.1, dynamicsymbols('x'): 0.0}


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
        # Modifying hte dict directly does change the dict.
        sys = System(self.kane)
        sys.initial_conditions[dynamicsymbols('x')] = 5.1
        assert sys.initial_conditions == dict(
                zip(dynamicsymbols('x, v'), [5.1, 1.0]))
        # Putting in a non-state key does not raise exception.
        sys.initial_conditions[symbols('m')] = 3.5
        # Then, if we integrate, we do error-checking and we get an exception.
        testing.assert_raises(ValueError, sys.integrate, [0.0, 1.0])


        # Keys must be coords or speeds.
        # ------------------------------
        testing.assert_raises(ValueError,
                self.sys.initial_conditions, {symbols('k'): 0.4}
                )
        testing.assert_raises(ValueError,
                self.sys.initial_conditions, {symbols('F'): 7.3}
                )

    def test_generate_ode_function(self):

        rhs = self.sys.generate_ode_function()

        assert rhs is self.sys.evaluate_ode_function

        args = {'constants': self.sys.constants,
                'specified': self.sys.specifieds}

        actual = rhs(np.ones(2), 0.0, args=(args,))

        print actual
        raise Exception()
        # TODO testing.assert_allclose(actual, np.array([1, 1]))
        # make a regression test.


        # n-link cart: play with specifieds.
        # ----------------------------------
        km = n_link_cart(3, cart_force=True, joint_torques=True, only_return_kane=True)
        sys = System(km)
        spec_syms = sys.specifieds.keys()
        rhs = sys.generate_ode_function()
        x = np.array(np.random.random(len(self.sys.states)))
        args = {'constants': {k: 1.0 for k in constants}}

        # Specify constants in two different ways and ensure we get the
        # same results. This is like Jason's test in codegen.
        args['specified'] = dict(zip(spec_syms, [1.0, 2.0, 3.0, 4.0]))
        xd_01 = rhs(x, 0.0, args)

        args['specified'] = {spec_syms[0]: lambda x, t: np.ones(t),
                (spec_syms[3], spec_syms[1]): lambda x, t: np.array([4, 2]),
                spec_syms[2]: 3.0 * np.ones(1)}
        xd_02 = rhs(x, 0.0, args)

        testing.assert_allclose(xd_01, xd_02)


    def test_integrate(self):

        times = np.linspace(0, 1, 100)

        # Try without calling generate_ode_function.
        # ------------------------------------------
        x_01 = self.sys.integrate(times)

        sys = System(self.kane)
        rhs = sys.generate_ode_function(generator='lambdify')
        x_02 = sys.integrate(times)

        testing.assert_allclose(x_01, x_02)

        # Ensure that initial conditions are reordered properly.
        # ------------------------------------------------------
        # I know that this is the order of the states.
        x0 = [5.1, 3.7]
        ic = {dynamicsymbols('x'): x[0], dynamicsymbols('v'): x[1]}
        sys.initial_conditions = ic
        x_03 = sys.integrate(times)
        x_04 = sys.ode_solver(sys.evaluate_ode_function,
                x0, times, args={'constants': sys.constants,
                    'specified': sys.specifieds})

        testing.assert_allclose(x_03, x_04)

        # Test a generator other than lambdify.
        # -------------------------------------
        # TODO

        # Unrecognized generator.
        # -----------------------
        sys = System(self.kane)
        sys.generate_ode_function(generator='made-up')


        raise Exception()








