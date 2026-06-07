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
import pytest
Cython = sm.external.import_module('Cython')
symjit = sm.external.import_module('symjit')

from ..system import System
from ..models import multi_mass_spring_damper, n_link_pendulum_on_cart
from ..utils import PyDyImportWarning, _sort_velocity_constraints

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

        assert sys.auxiliaries == []
        assert sys.constants == dict()
        assert sys.constants_symbols == set(sm.symbols('k0, m0, g, c0'))
        assert sys.constraints == sm.Matrix([])
        assert sys.coordinates == [me.dynamicsymbols('x0')]
        assert sys.eom_method is self.kane
        assert sys.evaluate_ode_function is None
        assert sys.initial_conditions == dict()
        assert sys.noncontributing_forces == []
        assert sys.num_auxiliaries == 0
        assert sys.num_config_constraints == 0
        assert sys.num_constants == 4
        assert sys.num_constraints == 0
        assert sys.num_coordinates == 1
        assert sys.num_motion_constraints == 0
        assert sys.num_outputs == 0
        assert sys.num_specifieds == 1
        assert sys.num_speeds == 1
        assert sys.num_states == 2
        assert sys.ode_solver is odeint
        assert sys.specifieds == dict()
        assert sys.specifieds_symbols == {self.specified_symbol}
        assert sys.speeds == [me.dynamicsymbols('v0')]
        assert sys.states == list(me.dynamicsymbols('x0, v0'))
        np.testing.assert_allclose(sys.times, np.array([]))

        # Set constants and specified symbols manually
        # --------------------------------------------
        sys = System(self.kane,
                     constants_symbols=sm.symbols('k0, m0, g, c0'),
                     specifieds_symbols=(self.specified_symbol,))
        assert sys.constants_symbols == set(sm.symbols('k0, m0, g, c0'))
        assert sys.specifieds_symbols == {self.specified_symbol}

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
        # TODO : This may not work with outputs.
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

        with pytest.raises(ValueError):
            self.sys.set_dependent_initial_conditions()

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

    def test_evaluate_ode(self):

        self.sys.times = np.array([0.0, 1.0])
        x = self.sys.evaluate_ode()
        x_expected = np.array([0.0, 10.3])
        np.testing.assert_allclose(x, x_expected)

        def force(x, t):
            f = 5.8*t
            return f

        # make sure t is passed through
        self.sys.specifieds = {self.specified_symbol: force}
        self.sys.initial_conditions = {self.sys.states[0]: 5.1,
                                       self.sys.states[1]: -4.5}
        self.sys.times = np.array([1.2, 1.3])
        x = self.sys.evaluate_ode()
        x_expected = np.array([-4.5, 10.58])
        np.testing.assert_allclose(x, x_expected)

        xdot = self.sys.evaluate_ode(x=[1.0, 2.0])
        assert xdot.shape == (2,)
        xdot = self.sys.evaluate_ode(x=[1.0, 2.0], t=3.0)
        assert xdot.shape == (2,)
        xdot = self.sys.evaluate_ode(x=[[0.2, 1.0], [0.2, 3.0], [0.2, 0.3]],
                                     t=[3.0, 5.0, 6.0])
        assert xdot.shape == (3, 2)
        with pytest.raises(ValueError):
            xdot = self.sys.evaluate_ode(x=[1.0, 2.0, 3.0])
        with pytest.raises(ValueError):
            xdot = self.sys.evaluate_ode(x=[1.0, 2.0], t=[1.0, 2.0])
        with pytest.raises(ValueError):
            xdot = self.sys.evaluate_ode(x=[[0.2, 1.0],
                                            [0.2, 3.0],
                                            [0.2, 0.3]],
                                         t=[3.0, 5.0])

    def test_evaluate_constraints(self):
        with pytest.raises(ValueError):
            self.sys.evaluate_constraints()
        with pytest.raises(ValueError):
            self.sys.evaluate_config_constraints()
        with pytest.raises(ValueError):
            self.sys.evaluate_motion_constraints()

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
        if symjit:
            sys.generate_ode_function(generator='symjit')
            sys.times = times
            x_06 = sys.integrate()
            testing.assert_allclose(x_04, x_06)
        else:
            warnings.warn("Symjit was not found so the related tests are being"
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

    def test_c_contiguous(self):

        def ode_solver(f, x0, t, args=None):
            a = np.zeros((len(x0), len(x0)))
            for i in range(len(x0)):
                a[:, 0] = np.asarray(x0)
            #print('c', a[0, :].data.c_contiguous)
            #print('f', a[:, 0].data.f_contiguous)
            f(a[:, 0], t[0], *args)  # test with f contiguous array
            return odeint(f, a[:, 0], t, args=args)

        sys = System(self.kane_nlink, ode_solver=ode_solver)
        sys.times = np.linspace(0.0, 1.0)
        sys.generate_ode_function(generator='cython', force_c_contiguous=False)
        with pytest.raises(ValueError):
            sys.integrate()
        sys.generate_ode_function(generator='cython')
        sys.integrate()

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

    fr, frstar = kane.kanes_equations(bodies, loads)

    sys = System(kane)

    sys.specifieds = {(beta, beta.diff(), beta.diff().diff()):
                      lambda x, t: np.array([1.0, 1.0, 1.0])}

    sys.times = np.linspace(0, 10, 20)

    sys.integrate()


def test_system_with_constraints(plot=False):
    """Rolling disc. Start with disc suspended above a ground plane. Use a
    holonomic constraint to keep it rolling on the ground plane. Add no-slip
    rolling constraints for nonholonomic."""

    r, g, m = sm.symbols('r, g, m', real=True)

    x, y, z = me.dynamicsymbols('q1, q2, q3', real=True)
    yaw, roll, pitch = me.dynamicsymbols('q4, q5, q6', real=True)
    u1, u2, u3, u4, u5, u6 = me.dynamicsymbols('u1, u2, u3, u4, u5, u6',
                                               real=True)

    N, A, B, C = sm.symbols('N, A, B, C', cls=me.ReferenceFrame)
    O, P, Q = sm.symbols('O, P, Q', cls=me.Point)

    A.orient_axis(N, yaw, N.z)
    B.orient_axis(A, roll, A.x)
    C.orient_axis(B, pitch, B.y)
    N_w_C = C.ang_vel_in(N)

    kd_eqs=(
        x.diff() - u1,
        y.diff() - u2,
        z.diff() - u3,
        N_w_C.dot(C.z) - u4,
        N_w_C.dot(C.x) - u5,
        N_w_C.dot(C.y) - u6,
    )

    A.set_ang_vel(N, u4*N.z)  # yaw frame
    B.set_ang_vel(A, u5*A.x)  # roll frame
    C.set_ang_vel(B, u6*B.y)  # pitch frame (the disc)

    P.set_pos(O, x*N.x + y*N.y + z*N.z)  # bottom of disc
    Q.set_pos(P, r*B.z)  # disc center

    O.set_vel(N, 0)
    P.set_vel(N, u1*N.x + u2*N.y + u3*N.z)
    Q.v2pt_theory(P, N, B)

    # holonomic constraint to force bottom point to be in xy plane
    holonomic = (Q.pos_from(O).dot(N.z) - r*sm.cos(roll),)

    # velocity of point fixed in disc at ground, should have velocity of zero
    v = Q.vel(N) + C.ang_vel_in(N).cross(-r*B.z)
    nonholonomic = (v.dot(A.x), v.dot(A.y), v.dot(N.z))

    inertia = (me.inertia(C, m*r**2/4, m*r**2/2, m*r**2/4), Q)
    disc = me.RigidBody('disc', Q, C, m, inertia)

    gravity = (Q, -m*g*N.z)

    q = [x, y, yaw, roll, pitch]
    qdot_exprs = list(sm.solve(kd_eqs, [x.diff(), y.diff(), yaw.diff(), roll.diff(),
                                   pitch.diff()]).values())
    h_idxs, nh_idxs = _sort_velocity_constraints(nonholonomic, q, qdot_exprs)
    holonomic_constraints = nonholonomic[h_idxs]
    nonholonomic_constraints = nonholonomic[nh_idxs]

    kane = me.KanesMethod(
        N,
        (x, y, yaw, roll, pitch),
        (u4, u5, u6),
        kd_eqs=kd_eqs,
        q_dependent=(z,),
        u_dependent=(u1, u2, u3),
        configuration_constraints=holonomic,
        velocity_constraints=nonholonomic,
        bodies=(disc,),
        forcelist=(gravity,),
    )
    fr, frstar = kane.kanes_equations()

    sys = System(kane)

    assert sys.num_config_constraints == 1
    assert sys.num_motion_constraints == 3
    assert sys.num_constraints == 4

    sys.constants = {
        r: 0.3,
        g: 9.81,
        m: 1.25,
    }

    speed = 10.0
    yaw0 = np.deg2rad(10.0)
    roll0 = np.deg2rad(10.0)

    z_guess, u1_guess, u2_guess, u3_guess = 0.1, 10.0, 1.0, 0.1

    sys.initial_conditions = {
        x: 1.0,
        y: -1.0,
        z: z_guess,
        yaw: yaw0,
        roll: roll0,
        pitch: 0.0,
        u1: u1_guess,  #speed*np.cos(yaw0),
        u2: u2_guess,  #speed*np.sin(yaw0),
        u3: u3_guess,
        u4: 0.0,
        u5: np.deg2rad(100.0),
        u6: speed/sys.constants[r],
    }

    xdot0 = sys.evaluate_ode()  # assumes t = 0.0

    sys.times = np.array([1.0, 2.0])

    sys.set_dependent_initial_conditions()
    x0 = sys.initial_conditions
    np.testing.assert_allclose(x0[x], 1.0)
    np.testing.assert_allclose(x0[y], -1.0)
    np.testing.assert_allclose(x0[z], 0.0, atol=1e-12)
    np.testing.assert_allclose(x0[yaw], yaw0)
    np.testing.assert_allclose(x0[roll], roll0)
    np.testing.assert_allclose(x0[pitch], 0.0, atol=1e-12)
    np.testing.assert_allclose(x0[u1], speed*np.cos(yaw0))
    np.testing.assert_allclose(x0[u2], speed*np.sin(yaw0))
    np.testing.assert_allclose(x0[u3], 0.0, atol=1e-12)
    np.testing.assert_allclose(x0[u4], 0.0, atol=1e-12)
    np.testing.assert_allclose(x0[u5], np.deg2rad(100.0))
    np.testing.assert_allclose(x0[u6], speed/sys.constants[r])

    np.testing.assert_allclose(sys.evaluate_config_constraints(), 0.0,
                               atol=1e-12)
    np.testing.assert_allclose(sys.evaluate_motion_constraints(), 0.0,
                               atol=1e-12)
    np.testing.assert_allclose(sys.evaluate_constraints(), 0.0, atol=1e-12)

    xdot0 = sys.evaluate_ode()
    xdot0_expected = np.array([
        9.84807753,
        1.7364817766693035,
        0.0,
        1.7453292519943295,
        33.33333333,
        0.0,
        -118.15025127,
        4.542636327766898,
        20.51657582308824,
        6.0614648807520375,
        1.0687998010934021,
        0.0,
    ])
    np.testing.assert_allclose(xdot0, xdot0_expected, rtol=1e-10, atol=1e-10)

    xdot0 = sys.evaluate_ode(x=np.ones(len(sys.states)))
    xdot0_expected = np.array([
        1.0,
        1.0,
        -0.5574077246549022,
        1.3817732906760363,
        1.4690424270048894,
        1.0,
        -3.7016314353618576,
        21.57392570087562,
        2.4498490534935637,
        0.9545054524443378,
        0.06103534404727443,
        0.0,
    ])
    np.testing.assert_allclose(xdot0, xdot0_expected, rtol=1e-10, atol=1e-10)

    # nonholonomic function of: {u1(t), u2(t), u3(t), u4(t), u6(t)}
    # u3 = 0
    z_guess, u2_guess, u3_guess, u6_guess = 0.0, 0.0, 0.0, 30.0
    sys.initial_conditions = {
        x: 1.0,
        y: -1.0,
        z: z_guess,
        yaw: yaw0,
        roll: roll0,
        pitch: 0.0,
        u1: speed*np.cos(yaw0),
        u2: u2_guess,  #speed*np.sin(yaw0),
        u3: u3_guess,  # 0.0
        u4: 0.0,
        u5: np.deg2rad(100.0),
        u6: u6_guess,  #speed/sys.constants[r],
    }

    with pytest.raises(ValueError):  # too many dep vars
        sys.set_dependent_initial_conditions(dep_vars=(z, u2, u3, u6, u4))

    sys.set_dependent_initial_conditions(dep_vars=(z, u2, u3, u6),
                                         use_jac=True, tol=1e-11)

    x0 = sys.initial_conditions
    np.testing.assert_allclose(x0[x], 1.0)
    np.testing.assert_allclose(x0[y], -1.0)
    np.testing.assert_allclose(x0[z], 0.0, atol=1e-10)
    np.testing.assert_allclose(x0[yaw], yaw0)
    np.testing.assert_allclose(x0[roll], roll0)
    np.testing.assert_allclose(x0[pitch], 0.0, atol=1e-10)
    np.testing.assert_allclose(x0[u1], speed*np.cos(yaw0))
    np.testing.assert_allclose(x0[u2], speed*np.sin(yaw0))
    np.testing.assert_allclose(x0[u3], 0.0, atol=1e-10)
    np.testing.assert_allclose(x0[u4], 0.0, atol=1e-10)
    np.testing.assert_allclose(x0[u5], np.deg2rad(100.0))
    np.testing.assert_allclose(x0[u6], speed/sys.constants[r])

    np.testing.assert_allclose(sys.evaluate_config_constraints(), 0.0,
                               atol=1e-10)
    np.testing.assert_allclose(sys.evaluate_motion_constraints(), 0.0,
                               atol=1e-10)
    np.testing.assert_allclose(sys.evaluate_constraints(), 0.0, atol=1e-10)

    fps = 30  # frames per second
    duration = 24.0  # seconds
    sys.times = np.linspace(0.0, duration, num=int(duration*fps))

    trajectories = sys.integrate()
    con_traj = sys.evaluate_constraints(x=trajectories)

    if plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(len(sys.states), 1, sharex=True,
                                 layout='constrained')
        for ax, traj, s in zip(axes, trajectories.T, sys.states):
            ax.plot(sys.times, traj)
            ax.set_ylabel(s)

        fig, ax = plt.subplots()
        ax.plot(trajectories[:, 0], trajectories[:, 1])
        ax.set_aspect('equal')

        fig, axes = plt.subplots(con_traj.shape[1], 1, sharex=True,
                                 layout='constrained')
        for ax, traj in zip(axes, con_traj.T):
            ax.plot(sys.times, traj)

        plt.show()


def test_system_with_noncontributing_forces(plot=False):

    ###########################################################################
    # Define the symbolic dynamics: double simple pendulum with damping and
    # noncontributing tension forces exposed, based on:
    # https://moorepants.github.io/learn-multibody-dynamics/noncontributing.html
    ###########################################################################
    m1, m2, l1, l2, c, g = sm.symbols('m1, m2, l1, l2, c, g')
    q1, q2, u1, u2 = me.dynamicsymbols('q1, q2, u1, u2')
    u3, u4, T1, T2 = me.dynamicsymbols('u3, u4, T1, T2')
    # TODO : Should these be dynamicsymbols to be consistent?
    int_T1, int_T2 = sm.Symbol('∫ T1 dt'), sm.Symbol('∫ T2 dt')

    N, A, B = sm.symbols('N, A, B', cls=me.ReferenceFrame)

    A.orient_axis(N, q1, N.z)
    B.orient_axis(N, q2, N.z)
    A.set_ang_vel(N, u1*N.z)
    B.set_ang_vel(N, u2*N.z)

    O_ = me.Point('O')
    P1 = O_.locatenew('P1', -l1*A.y)
    P2 = P1.locatenew('P2', -l2*B.y)

    O_.set_vel(N, 0)
    P1.v2pt_theory(O_, N, A)
    P1.set_vel(N, P1.vel(N) - u3*A.y)
    P2.v2pt_theory(P1, N, B)
    P2.set_vel(N, P2.vel(N) - u4*B.y)

    bob1 = me.Particle('bob1', P1, m1)
    bob2 = me.Particle('bob2', P2, m2)

    loads = (
        (P1, -m1*g*N.y + T1*A.y - T2*B.y),
        (P2, -m2*g*N.y + T2*B.y),
        (A, -c*u1*N.z + c*(u2 - u1)*N.z),
        (B, - c*(u2 - u1)*N.z),
    )

    kane = me.KanesMethod(
        N,
        (q1, q2),
        (u1, u2),
        kd_eqs=[q1.diff() - u1, q2.diff() - u2],
        bodies=(bob1, bob2),
        forcelist=loads,
        u_auxiliary=(u3, u4),
    )
    kane.kanes_equations()

    ###########################################################################
    # Check that the system will correctly build with automatic parsing of the
    # auxiliary equations.
    ###########################################################################
    sys = System(kane, noncontributing_forces=(T1, T2))
    assert sys._linear_idxs == [0, 1]
    assert sys._linear_outputs_symbols == [T1, T2]
    assert sys._num_linear_outputs == 2
    assert sys._num_simple_outputs == 0
    assert sys._simple_idxs == []
    assert sys._simple_outputs_symbols == []
    assert sys.auxiliaries == [int_T1, int_T2]
    assert sys.constants_symbols == {m1, m2, l1, l2, c, g}
    assert sys.noncontributing_forces == [T1, T2]
    assert sys.constraints == sm.Matrix([])
    assert sys.coordinates == [q1, q2]
    assert sys.num_constraints == 0
    assert sys.num_outputs == 2
    assert sys.outputs == {(T1, T2): kane.auxiliary_eqs}
    assert sys.outputs_symbols == [T1, T2]
    assert sys.speeds == [u1, u2]
    assert sys.states == [q1, q2, u1, u2, int_T1, int_T2]
    with pytest.raises(ValueError):
        sys.set_dependent_initial_conditions()
    np.testing.assert_allclose(sys.evaluate_outputs(), [2.0, 1.0])
    # TODO : If you haven't provided times it assumes time is always zero. Is
    # this a good idea?
    np.testing.assert_allclose(sys.evaluate_outputs(x=np.zeros((3, 6))),
                               np.array([[2.0, 1.0], [2.0, 1.0], [2.0, 1.0]]))
    np.testing.assert_allclose(sys.evaluate_outputs(x=np.zeros((3, 6)),
                                                    t=[1.0, 2.0, 3.0]),
                               np.array([[2.0, 1.0], [2.0, 1.0], [2.0, 1.0]]))
    with pytest.raises(ValueError):  # time should be a float
        np.testing.assert_allclose(sys.evaluate_outputs(t=[1.0, 2.0, 3.0]),
                                   np.array([[2.0, 1.0], [2.0, 1.0],
                                             [2.0, 1.0]]))
    with pytest.raises(ValueError):  # time should be an array
        np.testing.assert_allclose(sys.evaluate_outputs(x=np.zeros((3, 6)),
                                                        t=2.0),
                                   np.array([[2.0, 1.0], [2.0, 1.0],
                                             [2.0, 1.0]]))
    with pytest.raises(ValueError):  # x.shape[0] must equal len(t)
        np.testing.assert_allclose(sys.evaluate_outputs(x=np.zeros((3, 6)),
                                                        t=[1.0, 2.0]),
                                   np.array([[2.0, 1.0], [2.0, 1.0],
                                             [2.0, 1.0]]))

    ###########################################################################
    # Check a system with both constraints and noncontributing forces
    ###########################################################################
    # TODO : If u3 and u4 are not replaced here, equations of motion contain
    # derivatives of u3 and u4. Maybe KanesMethod should to this substitution
    # internally.
    T_ = me.dynamicsymbols('T')
    ke = (m1/2*P1.vel(N).dot(P1.vel(N)) +
          m2/2*P2.vel(N).dot(P2.vel(N))).xreplace({u3: 0, u4: 0})
    mot_con = P2.vel(N).dot(B.x).xreplace({u3: 0, u4: 0})
    kane_with_con = me.KanesMethod(
        N,
        (q1, q2),
        (u1,),
        kd_eqs=[q1.diff() - u1, q2.diff() - u2],
        u_dependent=(u2,),
        velocity_constraints=[mot_con],
        bodies=(bob1, bob2),
        forcelist=loads,
        u_auxiliary=(u3, u4),
    )
    kane_with_con.kanes_equations()

    ec = sm.symbols('ec')
    sys = System(
        kane_with_con,
        constants={m1: 1.0, m2: 2.0, l1: 1.0, l2: 2.0, c: 1.0, g: 9.81},
        initial_conditions={q1: 0.1, q2: -0.2, u1: 1.4},
        noncontributing_forces=(T1, T2),
        outputs={T_: ec*ke},  # make sure extra constant gets picked up
    )
    sys.set_dependent_initial_conditions()

    # TODO : Should the user have access to the dummy symbols for the
    # constraint residuals other than in sys.output_symbols?
    fn = sys._con_syms[0]
    assert sys.constants_symbols == {m1, m2, l1, l2, c, g, ec}
    assert sys._linear_outputs_symbols == [T1, T2]
    assert sys._num_linear_outputs == 2
    assert sys._num_simple_outputs == 2
    assert sys._simple_outputs_matrix == sm.Matrix([ec*ke, mot_con])
    assert sys._simple_outputs_symbols == [T_, fn]
    assert sys.noncontributing_forces == [T1, T2]
    assert sys.constraints == sm.Matrix([mot_con])
    assert sys.num_config_constraints == 0
    assert sys.num_constraints == 1
    assert sys.num_motion_constraints == 1
    assert sys.num_outputs == 4
    assert sys.outputs_symbols == [T_, fn, T1, T2]

    np.testing.assert_allclose(sys.evaluate_ode(),
        [1.4, -0.6687355423879242, -10.857684528809608, 5.614318317403088,
         30.235328069913542, 18.345324])

    np.testing.assert_allclose(sys.evaluate_outputs(),
        [1.151171,  0.0, 30.235328, 18.345324])

    np.testing.assert_allclose(sys.evaluate_motion_constraints(), [0.0])
    np.testing.assert_allclose(sys.evaluate_constraints(), [0.0])
    with pytest.raises(ValueError):
        sys.evaluate_config_constraints()

    ###########################################################################
    # Check that the system will skip automatic parsing of the auxiliary
    # equations and only handle additional simple outputs passed into __init__.
    ###########################################################################
    c_ = me.dynamicsymbols('c')
    constraint = P2.pos_from(O_).dot(N.x) - sm.sin(2)
    outputs = {T_: ke, c_: constraint}

    sys = System(kane, outputs=outputs)

    assert sys._linear_idxs == []
    assert sys._linear_outputs_symbols == []
    assert sys._num_linear_outputs == 0
    assert sys._num_simple_outputs == 2
    assert sys._simple_idxs == [0, 1]
    assert sys._simple_outputs_symbols == [T_, c_]
    assert sys.auxiliaries == []
    assert sys.noncontributing_forces == []
    assert sys.constraints == sm.Matrix([])
    assert sys.num_constraints == 0
    assert sys.num_outputs == 2
    assert sys.outputs == {T_: ke, c_: constraint}
    assert sys.outputs_symbols == [T_, c_]
    assert sys.states == [q1, q2, u1, u2]
    with pytest.raises(ValueError):
        sys.evaluate_constraints()
    with pytest.raises(ValueError):
        sys.evaluate_config_constraints()
    with pytest.raises(ValueError):
        sys.evaluate_motion_constraints()
    np.testing.assert_allclose(sys.evaluate_ode(),
                               [0.0, 0.0, 0.0, 0.0])
    np.testing.assert_allclose(sys.evaluate_outputs(),
                               [0.0, -0.9092974268256817])

    # cannot have duplicate output names
    with pytest.raises(ValueError):
        sys.outputs = {T_: ke, (T_, c_): (ke, constraint)}

    ###########################################################################
    # Now change the outputs to a new dictionary and make sure things update.
    ###########################################################################
    K_, c2, X1, X2, a = me.dynamicsymbols('K, c2, X1, X2, a')
    outputs = {
        T_: ke,
        X1: 2*X1 + 4*u1.diff() + 5*u2.diff() - 10,  # check singelton works
        c_: constraint,
        (K_, c2): [ke, constraint],
    }

    sys.outputs = outputs
    assert sys._needs_code_regeneration

    assert sys._linear_outputs_forcing_rows == sm.Matrix([10])
    assert sys._linear_outputs_mass_matrix_rows == sm.Matrix([[4, 5, 2]])
    assert sys._linear_outputs_symbols == [X1]
    assert sys._num_linear_outputs == 1

    outputs = {
        T_: ke,
        (X1, X2): sm.Matrix([2*X1 + 3*X2 + 4*u1.diff() + 5*u2.diff() - 10,
                             6*X1 + 7*X2 + 8*u1.diff() + 9*u2.diff() - 11]),
        c_: constraint,
        (K_, c2): sm.Matrix([ke, constraint]),
        a: a - u1.diff(),  # trick to return u'
    }

    sys.outputs = outputs
    assert sys._needs_code_regeneration

    assert sys._linear_outputs_forcing_rows == sm.Matrix([10, 11, 0])
    assert sys._linear_outputs_mass_matrix_rows == sm.Matrix([[4, 5, 2, 3, 0],
                                                              [8, 9, 6, 7, 0],
                                                              [-1, 0, 0, 0, 1]])
    assert sys._linear_outputs_symbols == [X1, X2, a]
    assert sys._num_linear_outputs == 3
    assert sys._num_simple_outputs == 4
    assert sys._simple_outputs_matrix == sm.Matrix([ke, constraint, ke,
                                                    constraint])
    assert sys._simple_outputs_symbols == [T_, c_, K_, c2]
    assert sys.num_outputs == 7
    assert sys.outputs_symbols == [T_, X1, X2, c_, K_, c2, a]

    np.testing.assert_allclose(sys.evaluate_outputs(),
        [0.0, -9.25, 9.5, -0.9092974268256817, 0.0, -0.9092974268256817, 0.0],
                               atol=1e-12)

    # Now when constraint loads are added, System will check KanesMethod for
    # any auxilliary equations for the noncontributing forces. These will be
    # appended to the outputs automatically. You have to add the correct number
    # of constraint load symbols and they have to be present in the auxiliary
    # equations.
    with pytest.raises(ValueError):  # too many symbols
        sys.noncontributing_forces = (T1, T2, u1)

    sys.noncontributing_forces = (T1, T2)
    assert sys._needs_code_regeneration

    assert sys._linear_outputs_symbols == [X1, X2, a, T1, T2]
    assert sys._num_linear_outputs == 5
    assert sys._num_simple_outputs == 4
    assert sys._simple_outputs_matrix == sm.Matrix([ke, constraint, ke,
                                                    constraint])
    assert sys.num_outputs == 9
    assert sys.outputs_symbols == [T_, X1, X2, c_, K_, c2, a, T1, T2]

    sys.constants = {
        m1: 1.0,
        m2: 2.0,
        l1: 1.0,
        l2: 2.0,
        c: 1.0,
        g: 9.81,
    }

    sys.initial_conditions = {
        q1: 0.1,
        q2: 0.0,
    }

    np.testing.assert_allclose(sys.evaluate_outputs(),
        [
            0.0,
            -9.26439137981748,
            10.96192493300433,
            -0.8094640101788535,
            0.0,
            -0.8094640101788535,
            -2.880676,
            28.710670665299798,
            19.044824599932618,
        ])

    sys.times = np.linspace(0.0, 4.0, num=400)

    x_traj = sys.integrate()
    assert x_traj.shape == (400, 9)  # shouldn't include dummy states?
    xdot_traj = sys.evaluate_ode(x=x_traj)
    assert x_traj.shape == (400, 9)  # shouldn't include dummy states?
    y_traj = sys.evaluate_outputs(x=x_traj)
    assert x_traj.shape == (400, 9)

    rhs = sys.generate_ode_function()
    print(rhs.__doc__)

    if plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(len(sys.states), 1, sharex=True,
                                 layout='constrained')
        for ax, traj, s in zip(axes, x_traj.T, sys.states):
            ax.plot(sys.times, traj)
            ax.set_ylabel(sm.latex(s, mode='inline'))

        fig, axes = plt.subplots(len(sys.states), 1, sharex=True,
                                 layout='constrained')
        for ax, traj, s in zip(axes, xdot_traj.T,
                               [xi.diff() for xi in sys.coordinates +
                                sys.speeds] +
                               sys._linear_outputs_symbols):
            ax.plot(sys.times, traj)
            ax.set_ylabel(sm.latex(s, mode='inline'))

        fig, axes = plt.subplots(sys.num_outputs, 1, sharex=True,
                                 layout='constrained')
        for ax, traj, s in zip(axes, y_traj.T, sys.outputs_symbols):
            ax.plot(sys.times, traj)
            ax.set_ylabel(sm.latex(s, mode='inline'))

        plt.show()


def test_explicit_time():
    sys = multi_mass_spring_damper(1, apply_gravity=True,
                                   apply_external_forces=True,
                                   external_force_as_sin_of_t=True)
    constant_map = dict(zip(sm.symbols('m0, k0, c0, g'),
                            [2.0, 1.5, 0.5, 9.8]))
    sys.constants = constant_map
    ts = np.linspace(0.0, 2.0, num=3)
    xd = sys.evaluate_ode()
    exp_xd = np.array([0.0, (-0.5*0.0 - 1.5*0.0 + 2.0*9.8)/2.0 +
                       np.sin(ts[0])/2.0])
    np.testing.assert_allclose(xd, exp_xd)
    xds = sys.evaluate_ode(x=np.ones((3, 2)), t=ts)
    exp_xds = np.array([[1.0, (-0.5*1.0 - 1.5*1.0 + 2.0*9.8)/2.0 +
                         np.sin(ts[0])/2.0],
                        [1.0, (-0.5*1.0 - 1.5*1.0 + 2.0*9.8)/2.0 +
                         np.sin(ts[1])/2.0],
                        [1.0, (-0.5*1.0 - 1.5*1.0 + 2.0*9.8)/2.0 +
                         np.sin(ts[2])/2.0]])
    np.testing.assert_allclose(xds, exp_xds)
    sys.times = ts
    xs = sys.integrate()
    exp_xs = np.array([[ 0.0, 0.0],
                       [ 4.31701656, 7.82602082],
                       [13.45662943, 9.24246248]])
    np.testing.assert_allclose(xs, exp_xs)


def test_issue_427():
    #  ensure integrate works with constants being an empty dict

    f = 3.0

    t = sm.Symbol('t')
    x, xv = me.dynamicsymbols('x xv')
    thetav = me.dynamicsymbols('thetav')
    xd = me.dynamicsymbols('x', 1)

    BaseFrame = me.ReferenceFrame('BaseFrame')
    RodFrame = BaseFrame.orientnew('RodFrame', 'Axis', [0, BaseFrame.z])
    RodFrame.set_ang_vel(BaseFrame, thetav*BaseFrame.z)

    origin  = me.Point('origin')
    rodcentre  = origin.locatenew('rodcentre', x*BaseFrame.x)
    rodcentre.set_vel(BaseFrame, rodcentre.pos_from(origin).diff(t, BaseFrame))

    Izz = me.outer(RodFrame.z, RodFrame.z)*10
    rodbody = me.RigidBody(name="rod", masscenter=rodcentre, frame=RodFrame,
                           mass=10, inertia=(Izz, rodcentre))

    forces = [(rodcentre, f*RodFrame.x)]
    bodies = [rodbody]

    KM = me.KanesMethod(BaseFrame, q_ind=[x], u_ind=[xv], kd_eqs=[ xd - xv ])
    KM.kanes_equations(bodies, forces)

    sys = System(KM,
                constants={},
                initial_conditions={x:0., xv:10.},
                times=np.linspace(0.0, 6, 60))
    sys.integrate()
