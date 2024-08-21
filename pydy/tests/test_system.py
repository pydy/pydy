#!/usr/bin/env python

import os
import warnings
import tempfile
import shutil

import numpy as np
from numpy import testing
import sympy as sm
import sympy.physics.mechanics as me
from scipy.integrate import odeint
theano = sm.external.import_module('theano')
Cython = sm.external.import_module('Cython')

from ..system import System
from ..models import multi_mass_spring_damper, n_link_pendulum_on_cart
from ..utils import PyDyImportWarning, sympy_newer_than

SYMPY_VERSION = sm.__version__

warnings.simplefilter('once', PyDyImportWarning)


class TestSystem():

    def setup_method(self):

        # Create a simple system with one specified quantity.
        self.sys = multi_mass_spring_damper(1, apply_gravity=True,
                                            apply_external_forces=True)
        self.specified_symbol = next(iter(self.sys.specifieds_symbols))
        self.constant_map = dict(zip(sm.symbols('m0, k0, c0, g'),
                                     [2.0, 1.5, 0.5, 9.8]))
        self.sys.specifieds = {self.specified_symbol: np.ones(1)}
        self.sys.constants = self.constant_map
        self.kane = self.sys.eom_method

        # Create a system with multiple specified quantities.
        self.kane_nlink = n_link_pendulum_on_cart(3, cart_force=True,
                                                  joint_torques=True).eom_method

    def test_init(self):

        # Check defaults for most attributes.
        # -----------------------------------
        sys = System(self.kane)

        assert (sys.constants_symbols ==
                set(sm.symbols('k0, m0, g, c0')))
        assert sys.specifieds_symbols == {self.specified_symbol}
        assert sys.states == me.dynamicsymbols('x0, v0')
        assert sys.evaluate_ode_function is None
        assert sys.eom_method is self.kane
        assert sys.ode_solver is odeint
        assert sys.specifieds == dict()
        assert sys.initial_conditions == dict()
        assert sys.constants == dict()
        assert sys.times == list()

        # Specify a bunch of attributes during construction.
        # --------------------------------------------------
        ic = {me.dynamicsymbols('x0'): 3.6, me.dynamicsymbols('v0'): 4.3}
        sys = System(self.kane,
                     ode_solver=odeint,
                     specifieds={self.specified_symbol: np.ones(1)},
                     initial_conditions=ic,
                     constants=self.constant_map)

        assert sys.eom_method is self.kane
        assert list(sys.specifieds.keys()) == [me.dynamicsymbols('f0')]
        testing.assert_allclose(list(sys.specifieds.values()),
                                [np.ones(1)])
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(list(sys.initial_conditions.values()),
                                list(ic.values()))
        assert sys.constants.keys() == self.constant_map.keys()
        testing.assert_allclose(list(sys.constants.values()),
                                list(self.constant_map.values()))

        # Use old specifieds.
        # -------------------
        sys = System(self.kane,
                     ode_solver=odeint,
                     specifieds={'symbols': [self.specified_symbol],
                                 'values': np.ones(1)},
                     initial_conditions=ic,
                     constants=self.constant_map)

    def test_coordinates(self):
        assert self.sys.coordinates == self.kane.q[:]

    def test_speeds(self):
        assert self.sys.speeds == self.kane.u[:]

    def test_states(self):
        assert self.sys.states == self.kane.q[:] + self.kane.u[:]

    def test_constants(self):

        # User-specified numerical values for constants.
        constants = {sm.symbols('m0'): 3.0, sm.symbols('c0'): 4.5}

        # The System should fill in the remaining constants with the defaults.

        # Construct a new system with our constants at construction time.
        # ---------------------------------------------------------------
        sys = System(self.kane, constants=constants)

        assert sys.constants.keys() == constants.keys()
        testing.assert_allclose(list(sys.constants.values()),
                                list(constants.values()))

        # Set constants after construction.
        # ---------------------------------
        sys = System(self.kane)

        # All constants have the default value.
        assert sys.constants == dict()
        sys.constants = constants

        assert sys.constants.keys() == constants.keys()
        testing.assert_allclose(list(sys.constants.values()),
                                list(constants.values()))

        # Using the property as a dict.
        # -----------------------------
        sys = System(self.kane)
        # Modifying the dict directly does change the dict.
        sys.constants[sm.symbols('m0')] = 9.3
        assert list(sys.constants.keys()) == [sm.symbols('m0')]
        testing.assert_allclose(list(sys.constants.values()), [9.3])

        # Putting in a non-constant key does not raise exception.
        sys.constants[me.dynamicsymbols('v0')] = 9.8
        # Then, if we integrate, we do error-checking and we get an exception.
        sys.times = [0.0, 1.0]
        with testing.assert_raises(ValueError):
            sys.integrate()

        # Provide a constant that isn't actually a constant.
        # --------------------------------------------------
        with testing.assert_raises(ValueError):
            sys.constants = {me.dynamicsymbols('x0'): 1.3}
        with testing.assert_raises(ValueError):
            sys.constants = {me.dynamicsymbols('f0'): 1.8}

    def test_specifieds(self):

        sys = System(self.kane)
        assert sys.specifieds == dict()
        sys.specifieds = {me.dynamicsymbols('f0'): 5.9}
        assert list(sys.specifieds.keys()) == [me.dynamicsymbols('f0')]
        testing.assert_allclose(list(sys.specifieds.values()), [5.9])

        # Using the property as a dict.
        # -----------------------------
        # Modifying the dict directly does change the dict.
        sys.specifieds[me.dynamicsymbols('f0')] = 5.1
        assert list(sys.specifieds.keys()) == [me.dynamicsymbols('f0')]
        testing.assert_allclose(list(sys.specifieds.values()), [5.1])
        # Putting in a non-specified key does not raise exception.
        sys.specifieds[me.dynamicsymbols('v0')] = 3.5
        # Then, if we integrate, we do error-checking and we get an exception.
        sys.times = [0.0, 1.0]
        with testing.assert_raises(ValueError):
            sys.integrate()

        sys = System(self.kane)
        # Putting in a value of the wrong length does not raise exception.
        sys.specifieds[me.dynamicsymbols('f0')] = 3.1 * np.ones(2)
        # Then, if we integrate, we do error-checking and we get an exception.
        # TODO actually, this does not seem to throw an exception.
        # TODO with testing.assert_raises(ValueError):
        # TODO     sys.integrate([0.0, 1.0])

        # The specified symbol must exist in the equations of motion and not
        # be a state.
        # ------------------------------------------------------------------
        sys = System(self.kane)
        with testing.assert_raises(ValueError):
            sys.specifieds = {sm.symbols('m0'): 5.4}
        with testing.assert_raises(ValueError):
            sys.specifieds = {me.dynamicsymbols('x0'): 5.1}

        # Complex error-checking when using property as a dict.
        # -----------------------------------------------------
        sys = System(self.kane_nlink)
        spec_syms = list(sys.specifieds_symbols)
        times = np.linspace(0, 0.5, 10)
        sys.specifieds = {
            spec_syms[0]: lambda x, t: np.ones(t),
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
        # Also check here that optional arguments can be passed to the solver.
        sys.specifieds.pop(spec_syms[0])
        state_traj, infodict = sys.integrate(full_output=True)

        # Test old way of providing specifieds.
        # -------------------------------------
        sys = System(self.kane_nlink)
        spec_syms = list(sys.specifieds_symbols)
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
            sys.specifieds = {'symbols': [sm.symbols('m1')], 'values': [1.0]}
        with testing.assert_raises(ValueError):
            sys.specifieds = {'symbols': [sm.symbols('T2, T2')],
                              'values': [1, 2]}
        with testing.assert_raises(ValueError):
            sys.specifieds = {'symbols': [me.dynamicsymbols('T2')],
                              'values': [1.0]}

        # Reordering causes issues!
        # -------------------------
        sys.specifieds = {
            'symbols': [spec_syms[1], spec_syms[0], spec_syms[2],
                        spec_syms[3]],
            'values': [2.0, 1.0, 3.0, 4.0]}
        # I tested: x_01 is not allclose to x_03.

        sys.generate_ode_function()
        x_04 = sys.integrate()
        testing.assert_allclose(x_01, x_04)

        # Test with no specifieds.
        sys = multi_mass_spring_damper(1, apply_gravity=True)
        sys.initial_conditions = {me.dynamicsymbols('x0'): 0.1,
                                  me.dynamicsymbols('v0'): -1.0}
        sys.times = times
        sys.integrate()

    def test_times(self):
        times1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        times2 = [0, -2, 7, 3, -5]
        times3 = [1, 2, 7, 4, 5]
        times4 = 4

        sys = System(self.kane, times=times1)
        testing.assert_allclose(sys.times, times1)

        with testing.assert_raises(ValueError):
            sys.times = times2

        with testing.assert_raises(ValueError):
            sys.times = times3

        with testing.assert_raises(TypeError):
            sys.times = times4

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
        ic = {me.dynamicsymbols('v0'): 6.1}

        # Using the constructor.
        # ----------------------
        sys = System(self.kane, initial_conditions=ic)
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(list(sys.initial_conditions.values()),
                                list(ic.values()))

        # Set the attribute.
        # ------------------
        sys = System(self.kane)
        sys.initial_conditions = ic
        assert sys.initial_conditions.keys() == ic.keys()
        testing.assert_allclose(list(sys.initial_conditions.values()),
                                list(ic.values()))

        # Using the property as a dict.
        # -----------------------------
        # Modifying hte dict directly does change the dict.
        sys = System(self.kane, times=[0.0, 1.0])
        sys.initial_conditions[me.dynamicsymbols('x0')] = 5.8
        assert list(sys.initial_conditions.keys()) == [me.dynamicsymbols('x0')]
        testing.assert_allclose(list(sys.initial_conditions.values()), [5.8])
        # Putting in a non-state key does not raise exception.
        sys.initial_conditions[sm.symbols('m0')] = 7.9
        # Then, if we integrate, we do error-checking and we get an exception.

        with testing.assert_raises(ValueError):
            sys.integrate()

        # Keys must be coords or speeds.
        # ------------------------------
        with testing.assert_raises(ValueError):
                self.sys.initial_conditions = {sm.symbols('k0'): 0.4}
        with testing.assert_raises(ValueError):
                self.sys.initial_conditions = {sm.symbols('f0'): 7.3}

    def test_generate_ode_function(self):

        rhs = self.sys.generate_ode_function()

        assert rhs is self.sys.evaluate_ode_function

        args = (self.sys.specifieds, self.sys.constants)

        actual = rhs(np.ones(2), 0.0, *args)

        # Regression.
        testing.assert_allclose(actual, np.array([1, 9.3]))

        # n-link cart: play with specifieds.
        # ----------------------------------
        sys = System(self.kane_nlink)
        spec_syms = list(sys.specifieds_symbols)
        rhs = sys.generate_ode_function()
        x = np.array(np.random.random(len(sys.states)))
        args = (self.sys.specifieds,
                {k: 1.0 for k in sys.constants_symbols})

        # Specify constants in two different ways and ensure we get the
        # same results. This is like Jason's test in codegen.
        args = (dict(zip(spec_syms, [1.0, 2.0, 3.0, 4.0])),
                {k: 1.0 for k in sys.constants_symbols})
        xd_01 = rhs(x, 0.0, *args)

        args = ({spec_syms[0]: lambda x, t: np.ones(1),
                 (spec_syms[3], spec_syms[1]): lambda x, t: np.array([4, 2]),
                 spec_syms[2]: 3.0 * np.ones(1)},
                {k: 1.0 for k in sys.constants_symbols})
        xd_02 = rhs(x, 0.0, *args)

        testing.assert_allclose(xd_01, xd_02)

    def test_integrate(self):

        times = np.linspace(0, 1, 100)

        # Try without calling generate_ode_function.
        # ------------------------------------------
        sys = System(self.kane, times=times)
        x_01 = sys.integrate()

        sys = System(self.kane, times=times)
        sys.generate_ode_function(generator='lambdify')
        x_02 = sys.integrate()

        testing.assert_allclose(x_01, x_02)

        # Ensure that the defaults are as expected.
        # -----------------------------------------
        constants_dict = dict(zip(sm.symbols('m0, k0, c0, g'),
                                  [1.0, 1.0, 1.0, 1.0]))
        specified_dict = {me.dynamicsymbols('f0'): 0.0}
        x_03 = sys.ode_solver(sys.evaluate_ode_function, [0, 0], sys.times,
                              args=(specified_dict, constants_dict))
        testing.assert_allclose(x_02, x_03)

        # Ensure that initial conditions are reordered properly.
        # ------------------------------------------------------
        sys = System(self.kane, times=times)
        # I know that this is the order of the states.
        x0 = [5.1, 3.7]
        ic = {me.dynamicsymbols('x0'): x0[0], me.dynamicsymbols('v0'): x0[1]}
        sys.initial_conditions = ic
        x_04 = sys.integrate()
        x_05 = sys.ode_solver(
            sys.evaluate_ode_function, x0, sys.times,
            args=(sys._specifieds_padded_with_defaults(),
                  sys._constants_padded_with_defaults()))

        testing.assert_allclose(x_04, x_05)

        # Test a generator other than lambdify.
        # -------------------------------------
        if theano:
            sys.generate_ode_function(generator='theano')
            sys.times = times
            x_06 = sys.integrate()
            testing.assert_allclose(x_04, x_06)
        else:
            warnings.warn("Theano was not found so the related tests are being"
                          " skipped.", PyDyImportWarning)

        # Unrecognized generator.
        # -----------------------
        sys = System(self.kane, times=times)
        with testing.assert_raises(NotImplementedError):
            sys.generate_ode_function(generator='made-up')

        # Test pass kwargs to the generators.
        if Cython:
            self.tempdirpath = tempfile.mkdtemp()
            prefix = 'my_test_file'
            self.sys.generate_ode_function(generator='cython',
                                           prefix=prefix,
                                           tmp_dir=self.tempdirpath)
            assert [True for f in os.listdir(self.tempdirpath)
                    if f.startswith(prefix)]
        else:
            warnings.warn("Cython was not found so the related tests are being"
                          " skipped.", PyDyImportWarning)
    def cleanup(self):
        shutil.rmtree(self.tempdirpath)


def test_specifying_coordinate_issue_339():
    """This test ensures that you can use derivatives as specified values."""

    # beta will be a specified angle
    beta = me.dynamicsymbols('beta')
    q1, q2, q3, q4 = me.dynamicsymbols('q1, q2, q3, q4')
    u1, u2, u3, u4 = me.dynamicsymbols('u1, u2, u3, u4')

    N = me.ReferenceFrame('N')
    A = N.orientnew('A', 'Axis', (q1, N.x))
    B = A.orientnew('B', 'Axis', (beta, A.y))

    No = me.Point('No')
    Ao = No.locatenew('Ao', q2 * N.x + q3 * N.y + q4 * N.z)
    Bo = Ao.locatenew('Bo', 10 * A.x + 10 * A.y + 10 * A.z)

    A.set_ang_vel(N, u1 * N.x)
    B.ang_vel_in(N)  # compute it automatically

    No.set_vel(N, 0)
    Ao.set_vel(N, u2 * N.x + u3 * N.y + u4 * N.z)
    Bo.v2pt_theory(Ao, N, B)

    body_A = me.RigidBody('A', Ao, A, 1.0, (me.inertia(A, 1, 2, 3), Ao))
    body_B = me.RigidBody('B', Bo, B, 1.0, (me.inertia(A, 3, 2, 1), Bo))

    bodies = [body_A, body_B]
    # TODO : This should be able to be simple an empty iterable.
    loads = [(No, 0 * N.x)]

    kdes = [u1 - q1.diff(),
            u2 - q2.diff(),
            u3 - q3.diff(),
            u4 - q4.diff()]

    kane = me.KanesMethod(N, q_ind=[q1, q2, q3, q4],
                          u_ind=[u1, u2, u3, u4], kd_eqs=kdes)

    if sympy_newer_than('1.0'):
        fr, frstar = kane.kanes_equations(bodies, loads)
    else:
        fr, frstar = kane.kanes_equations(loads, bodies)

    sys = System(kane)

    sys.specifieds = {(beta, beta.diff(), beta.diff().diff()):
                      lambda x, t: np.array([1.0, 1.0, 1.0])}

    sys.times = np.linspace(0, 10, 20)

    sys.integrate()
