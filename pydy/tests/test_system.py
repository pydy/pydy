
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

        self.kane_nlink = n_link_cart(3, cart_force=True, joint_torques=True,
                only_return_kane=True)

    def test_init(self):

        # Check defaults for most attributes.
        # -----------------------------------
        sys = System(self.kane)

        assert sys.constants_symbols == list(symbols('k, c, m, g'))
        assert sys.specifieds_symbols == [self.specified_symbol]
        assert sys.states == dynamicsymbols('x, v')
        assert sys.evaluate_ode_function == None
        assert sys.eom_method is self.kane
        assert sys.ode_solver is odeint
        assert sys.specifieds == dict()
        assert sys.initial_conditions == dict()
        assert sys.constants == dict()

        # Specify a bunch of attributes during construction.
        # --------------------------------------------------
        ic = {dynamicsymbols('x'): 3.6, dynamicsymbols('v'): 4.3}
        sys = System(self.kane,
                     ode_solver=odeint,
                     specifieds={self.specified_symbol: np.ones(1)},
                     initial_conditions=ic,
                     constants=self.constant_map)

        assert sys.eom_method is self.kane
        assert sys.specifieds.keys() == [dynamicsymbols('F')]
        testing.assert_allclose(sys.specifieds.values(), [np.ones(1)])
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(sys.initial_conditions.values(), ic.values())
        assert sys.constants.keys() == self.constant_map.keys()
        testing.assert_allclose(sys.constants.values(),
                self.constant_map.values())

        # Use old specifieds.
        # -------------------
        sys = System(self.kane,
                     ode_solver=odeint,
                     specifieds={'symbols': [self.specified_symbol],
                         'values': np.ones(1)},
                     initial_conditions=ic,
                     constants=self.constant_map)

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

        # Construct a new system with our constants at construction time.
        # ---------------------------------------------------------------
        sys = System(self.kane, constants=constants)

        assert sys.constants.keys() == constants.keys()
        testing.assert_allclose(sys.constants.values(), constants.values())


        # Set constants after construction.
        # ---------------------------------
        sys = System(self.kane)

        # All constants have the default value.
        assert sys.constants == dict()
        sys.constants = constants

        assert sys.constants.keys() == constants.keys()
        testing.assert_allclose(sys.constants.values(), constants.values())


        # Using the property as a dict.
        # -----------------------------
        sys = System(self.kane)
        # Modifying the dict directly does change the dict.
        sys.constants[symbols('m')] = 9.3
        assert sys.constants.keys() == [symbols('m')]
        testing.assert_allclose(sys.constants.values(), [9.3])

        # Putting in a non-constant key does not raise exception.
        sys.constants[dynamicsymbols('v')] = 9.8
        # Then, if we integrate, we do error-checking and we get an exception.
        sys.times = [0.0, 1.0]
        with testing.assert_raises(ValueError):
            sys.integrate()


        # Provide a constant that isn't actually a constant.
        # --------------------------------------------------
        with testing.assert_raises(ValueError):
            sys.constants = {dynamicsymbols('x'): 1.3}
        with testing.assert_raises(ValueError):
            sys.constants = {dynamicsymbols('F'): 1.8}

    def test_specifieds(self):

        sys = System(self.kane)
        assert sys.specifieds == dict()
        sys.specifieds = {dynamicsymbols('F'): 5.9}
        assert sys.specifieds.keys() == [dynamicsymbols('F')]
        testing.assert_allclose(sys.specifieds.values(), [5.9])

        # Using the property as a dict.
        # -----------------------------
        # Modifying the dict directly does change the dict.
        sys.specifieds[dynamicsymbols('F')] = 5.1
        assert sys.specifieds.keys() == [dynamicsymbols('F')]
        testing.assert_allclose(sys.specifieds.values(), [5.1])
        # Putting in a non-specified key does not raise exception.
        sys.specifieds[dynamicsymbols('v')] = 3.5
        # Then, if we integrate, we do error-checking and we get an exception.
        sys.times = [0.0, 1.0]
        with testing.assert_raises(ValueError):
            sys.integrate()

        sys = System(self.kane)
        # Putting in a value of the wrong length does not raise exception.
        sys.specifieds[dynamicsymbols('F')] = 3.1 * np.ones(2)
        # Then, if we integrate, we do error-checking and we get an exception.
        # TODO actually, this does not seem to throw an exception.
        # TODO with testing.assert_raises(ValueError):
        # TODO     sys.integrate([0.0, 1.0])


        # The specified symbol must exist in the equations of motion and not
        # be a state.
        # ------------------------------------------------------------------
        sys = System(self.kane)
        with testing.assert_raises(ValueError):
            sys.specifieds = {symbols('m'): 5.4}
        with testing.assert_raises(ValueError):
            sys.specifieds = {dynamicsymbols('x'): 5.1}

        # Complex error-checking when using property as a dict.
        # -----------------------------------------------------
        sys = System(self.kane_nlink)
        spec_syms = sys.specifieds_symbols
        times = np.linspace(0, 0.5, 10)
        sys.specifieds = {spec_syms[0]: lambda x, t: np.ones(t),
                (spec_syms[3], spec_syms[1]): lambda x, t: np.array([4, 2]),
                spec_syms[2]: 3.0 * np.ones(1)}
        # These won't throw an exception b/c we're modifying the dict directly.
        sys.specifieds[spec_syms[1]] = 7.1
        sys.times = times
        # This does.
        with testing.assert_raises(ValueError):
            sys.integrate()

        sys = System(self.kane_nlink)
        # This puts too many entries in the dict.
        sys.specifieds[spec_syms[0]] = 3.7
        sys.specifieds[(spec_syms[0], spec_syms[3])] = 5.8
        sys.times = times
        with testing.assert_raises(ValueError):
            sys.integrate()

        # This gets rid of the previous default entries, and should work
        # properly.
        sys.specifieds.pop(spec_syms[0])
        sys.integrate()

        # Test old way of providing specifieds.
        # -------------------------------------
        sys = System(self.kane_nlink)
        spec_syms = sys.specifieds_symbols
        # Get numbers using the new way.
        sys.specifieds = dict(zip(spec_syms, [1.0, 2.0, 3.0, 4.0]))
        sys.times = times
        x_01 = sys.integrate()

        # Now use the old way.
        sys.specifieds = {'symbols': spec_syms,
                'values': [1.0, 2.0, 3.0, 4.0]}
        x_02 = sys.integrate()
        testing.assert_allclose(x_01, x_02)

        # Error checks for the new way.
        # -----------------------------
        with testing.assert_raises(ValueError):
            sys.specifieds = {'symbols': [symbols('m1')], 'values': [1.0]}
        with testing.assert_raises(ValueError):
            sys.specifieds = {'symbols': [symbols('T2, T2')], 'values': [1, 2]}
        with testing.assert_raises(ValueError):
            sys.specifieds = {'symbols': [dynamicsymbols('T2')], 'values':
                    [1.0]}

        # Reordering causes issues!
        # -------------------------
        sys.specifieds = {'symbols':
                [spec_syms[1], spec_syms[0], spec_syms[2], spec_syms[3]],
                'values': [2.0, 1.0, 3.0, 4.0]}
        x_03 = sys.integrate()
        # I tested: x_01 is not allclose to x_03.

        sys.generate_ode_function()
        x_04 = sys.integrate()
        testing.assert_allclose(x_01, x_04)


    def test_ode_solver(self):

        assert self.sys.ode_solver == odeint
        self.sys.ode_solver = max
        assert self.sys.ode_solver is max

        # ode_solver must be a function
        # -----------------------------
        with testing.assert_raises(ValueError):
            self.sys.ode_solver = 5

    def test_initial_conditions(self):

        # Partially provided ic's.
        ic = {dynamicsymbols('v'): 6.1}


        # Using the constructor.
        # ----------------------
        sys = System(self.kane, initial_conditions=ic)
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(sys.initial_conditions.values(), ic.values())


        # Set the attribute.
        # ------------------
        sys = System(self.kane)
        sys.initial_conditions = ic
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(sys.initial_conditions.values(), ic.values())


        # Using the property as a dict.
        # -----------------------------
        # Modifying hte dict directly does change the dict.
        sys = System(self.kane, times=[0.0, 1.0])
        sys.initial_conditions[dynamicsymbols('x')] = 5.8
        assert sys.initial_conditions.keys() == [dynamicsymbols('x')]
        testing.assert_allclose(sys.initial_conditions.values(), [5.8])
        # Putting in a non-state key does not raise exception.
        sys.initial_conditions[symbols('m')] = 7.9
        # Then, if we integrate, we do error-checking and we get an exception.

        with testing.assert_raises(ValueError):
            sys.integrate()


        # Keys must be coords or speeds.
        # ------------------------------
        with testing.assert_raises(ValueError):
                self.sys.initial_conditions = {symbols('k'): 0.4}
        with testing.assert_raises(ValueError):
                self.sys.initial_conditions = {symbols('F'): 7.3}

    def test_generate_ode_function(self):

        rhs = self.sys.generate_ode_function()

        assert rhs is self.sys.evaluate_ode_function

        args = {'constants': self.sys.constants,
                'specified': self.sys.specifieds}

        actual = rhs(np.ones(2), 0.0, args)

        # Regression.
        testing.assert_allclose(actual, np.array([1, 9.3]))


        # n-link cart: play with specifieds.
        # ----------------------------------
        sys = System(self.kane_nlink)
        spec_syms = sys.specifieds_symbols
        rhs = sys.generate_ode_function()
        x = np.array(np.random.random(len(sys.states)))
        args = {'constants': {k: 1.0 for k in sys.constants_symbols}}

        # Specify constants in two different ways and ensure we get the
        # same results. This is like Jason's test in codegen.
        args['specified'] = dict(zip(spec_syms, [1.0, 2.0, 3.0, 4.0]))
        xd_01 = rhs(x, 0.0, args)

        args['specified'] = {spec_syms[0]: lambda x, t: np.ones(1),
                (spec_syms[3], spec_syms[1]): lambda x, t: np.array([4, 2]),
                spec_syms[2]: 3.0 * np.ones(1)}
        xd_02 = rhs(x, 0.0, args)

        testing.assert_allclose(xd_01, xd_02)


    def test_integrate(self):

        times = np.linspace(0, 1, 100)

        # Try without calling generate_ode_function.
        # ------------------------------------------
        sys = System(self.kane, times=times)
        x_01 = sys.integrate()

        sys = System(self.kane, times=times)
        rhs = sys.generate_ode_function(generator='lambdify')
        x_02 = sys.integrate()

        testing.assert_allclose(x_01, x_02)

        # Ensure that the defaults are as expected.
        # -----------------------------------------
        constants_dict = dict(zip(symbols('m, k, c, g'), [1.0, 1.0, 1.0, 1.0]))
        specified_dict = {dynamicsymbols('F'): 0.0}
        x_03 = sys.ode_solver(sys.evaluate_ode_function, [0, 0], sys.times,
                args=({
                    'constants': constants_dict,
                    'specified': specified_dict},))
        testing.assert_allclose(x_02, x_03)

        # Ensure that initial conditions are reordered properly.
        # ------------------------------------------------------
        sys = System(self.kane, times=times)
        # I know that this is the order of the states.
        x0 = [5.1, 3.7]
        ic = {dynamicsymbols('x'): x0[0], dynamicsymbols('v'): x0[1]}
        sys.initial_conditions = ic
        x_04 = sys.integrate()
        x_05 = sys.ode_solver(sys.evaluate_ode_function,
                x0, sys.times, args=(
                    {'constants': sys._constants_padded_with_defaults(),
                     'specified': sys._specifieds_padded_with_defaults()},))

        testing.assert_allclose(x_04, x_05)

        # Test a generator other than lambdify.
        # -------------------------------------
        sys.generate_ode_function(generator='theano')
        sys.times = times
        x_06 = sys.integrate()
        testing.assert_allclose(x_04, x_06)

        # Unrecognized generator.
        # -----------------------
        sys = System(self.kane, times=times)
        with testing.assert_raises(NotImplementedError):
            sys.generate_ode_function(generator='made-up')
