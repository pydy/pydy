"""The System class manages the simulation (integration) of a system whose
equations are given by
:external+sympy:py:class:`~sympy.physics.mechanics.kane.KanesMethod`.

Many of the attributes are also properties, and can be directly modified.

Here is the procedure for using this class.

    1. specify your options either via the constructor or via the
       attributes.
    2. optionally, call :py:meth:`~pydy.system.System.generate_ode_function` if
       you want to customize how the ODE function is generated.
    3. call :py:meth:`~pydy.system.System.integrate` to simulate your system.

The simplest usage of this class is as follows. First, we need a
:external+sympy:py:class:`~sympy.physics.mechanics.kane.KanesMethod` object on
which we have already invoked
:external+sympy:py:meth:`~sympy.physics.mechanics.kane.KanesMethod.kanes_equations`::

    >>> from sympy.physics.mechanics.models import n_link_pendulum_on_cart
    >>> from pydy.system import System
    >>> import numpy as np
    >>> kane = n_link_pendulum_on_cart()
    >>> times = np.linspace(0.0, 5.0, num=3)
    >>> sys = System(kane, times=times)
    >>> sys.integrate()
    array([[0., 0., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])

In this case, we use defaults for the numerical values of the constants,
specified quantities, initial conditions, etc. You probably won't like these
defaults. You can also specify such values via constructor keyword arguments or
via the attributes::

    >>> import sympy as sm
    >>> sys = System(kane,
    ...              initial_conditions={kane.q[1]: 0.5},
    ...              times=times)
    ...
    >>> g, l0, m0, m1 = list(sm.ordered(sys.constants_symbols))
    >>> sys.constants = {m1: 5.0}
    >>> sys.integrate()
    array([[ 0.        ,  0.5       ,  0.        ,  0.        ],
           [-1.12276473,  4.19253522, -0.77003647,  1.86016638],
           [-1.00443253,  5.47085374, -0.4536987 , -0.7915558 ]])

To double-check the constants, specifieds, states and times in your problem,
look at these properties::

    >>> sys.coordinates
    [q0(t), q1(t)]
    >>> sys.speeds
    [u0(t), u1(t)]
    >>> sys.states
    [q0(t), q1(t), u0(t), u1(t)]
    >>> sys.constants_symbols  # doctest: +SKIP
    {g, l0, m0, m1}
    >>> sys.specifieds_symbols
    {F(t)}
    >>> sys.times
    array([0. , 2.5, 5. ])

You can also add additional equations to evaluate alongside the differential
equations::

    >>> k0 = sm.Symbol('k0')
    >>> sys = System(kane,
    ...              initial_conditions={kane.q[1]: 0.5},
    ...              times=times,
    ...              outputs={k0: m0*sys.speeds[0]**2/2})
    >>> x = sys.integrate()
    >>> sys.outputs_symbols
    [k0]
    >>> sys.evaluate_outputs(x=x)
    array([[0.00000000e+00],
           [8.92619991e-01],
           [7.01894807e-04]])

In the prior examples, the :py:class:`System` generates the numerical ode
function for you behind the scenes. If you want to customize how this function
is generated, you must call
:py:meth:`~pydy.system.System.generate_ode_function` on your own::

    >>> rhs = sys.generate_ode_function(generator='cython')
    >>> sys.evaluate_ode_function == rhs
    True
    >>> help(sys.evaluate_ode_function)
    Help on function rhs in module pydy.codegen.ode_function_generators:
    <BLANKLINE>
    rhs(*args)
        Returns the derivatives of the states, i.e. numerically evaluates the right
        hand side of the first order differential equation.
    <BLANKLINE>
        x' = f(x, t, r, p)
    <BLANKLINE>
        Parameters
        ==========
        x : ndarray, shape(4,)
            The state vector is ordered as such:
                - q0(t)
                - q1(t)
                - u0(t)
                - u1(t)
        t : float
            The current time.
        r : dictionary; ndarray, shape(1,); function
    <BLANKLINE>
            There are three options for this argument. (1) is more flexible but
            (2) and (3) are much more efficient.
    <BLANKLINE>
            (1) A dictionary that maps the specified functions of time to floats,
            ndarrays, or functions that produce ndarrays. The keys can be a single
            specified symbolic function of time or a tuple of symbols. The total
            number of symbols must be equal to 1. If the value is a
            function it must be of the form g(x, t), where x is the current state
            vector ndarray and t is the current time float and it must return an
            ndarray of the correct shape. For example::
    <BLANKLINE>
              r = {a: 1.0,
                   (d, b) : np.array([1.0, 2.0]),
                   (e, f) : lambda x, t: np.array(x[0], x[1]),
                   c: lambda x, t: np.array(x[2])}
    <BLANKLINE>
            (2) A ndarray with the specified values in the correct order and of the
            correct shape.
    <BLANKLINE>
            (3) A function that must be of the form g(x, t), where x is the current
            state vector and t is the current time and it must return an ndarray of
            the correct shape.
    <BLANKLINE>
            The specified inputs are, in order:
                - F(t)
        p : dictionary len(4) or ndarray shape(4,)
            Either a dictionary that maps the constants symbols to their numerical
            values or an array with the constants in the following order:
                - g
                - l0
                - m0
                - m1
    <BLANKLINE>
        Returns
        =======
        dx : ndarray, shape(4,)
            The derivative of the state vector.
        y : ndarray, shape(1,)
            Values of the provided outputs.
                - y0(t)

    >>> sys.integrate()
    array([[ 0.        ,  0.5       ,  0.        ,  0.        ],
           [-0.31425675,  3.29123866, -1.33612873,  2.70246056],
           [-0.48148282,  5.77849021, -0.03746718, -0.08560791]])

"""
import warnings
from itertools import repeat

import numpy as np
import sympy as sm
from sympy.physics.mechanics import dynamicsymbols, find_dynamicsymbols
from scipy.integrate import odeint
from scipy.optimize import root

from .codegen.ode_function_generators import generate_ode_function
from .utils import (PyDyFutureWarning, PyDyUserWarning,
                    _sort_velocity_constraints)

SYMPY_VERSION = sm.__version__

warnings.simplefilter('once', PyDyFutureWarning)


class System(object):
    """Multibody dynamics system for simulation and numerical evaluation.

    See the class's attributes for a description of the arguments to this
    constructor.

    The parameters to this constructor are all attributes of the System. With
    the exception of :py:attr:`~pydy.system.System.eom_method`, these
    attributes can be modified directly at any future point.

    Parameters
    ----------
    eom_method : sympy.physics.mechanics.kane.KanesMethod
        You must have called
        :external+sympy:py:meth:`~sympy.physics.mechanics.kane.KanesMethod.kanes_equations`
        *before* constructing this system.
    constants : dict, optional (default: all 1.0)
        This dictionary maps SymPy
        :external+sympy:py:class:`~sympy.core.symbol.Symbol` objects to floats.
    specifieds : dict, optional (default: all 0.0)
        This dictionary maps SymPy Functions of time objects, or tuples of
        them, to floats, NumPy arrays, or functions of the state and time.
    ode_solver : function, optional
        This function computes the derivatives of the states. The default is
        :external+scipy:py:func:`scipy.integrate.odeint`.
    initial_conditions : dict, optional (default: all zero)
        This dictionary maps SymPy Functions of time objects to floats.
    times : array_like, shape(n,), optional
        An array_like object, which contains time values over which equations
        are integrated. It has to be supplied before
        :py:meth:`~System.integrate` can be called.
    outputs : dictionary, optional
        Maps functions of time or tuples of functions of time to expressions or
        iterables of expressions, respectively. In general, the expressions
        should be a function of the state, constants, and specifieds.
        Expressions that are linear in the functions of time and/or the time
        derivatives of the speeds are also supported, but not yet nonlinear
        functions of these variables.
    noncontributing_forces : iterable of Functions of time, optional
        If the ``eom_method`` includes noncontributig forces (Kane's method),
        provide a list of variable names for these forces and they will be
        computed when evaluating the differential equations.
    constants_symbols : iterable of Symbol, optional
        If provided, the system's equations will not be searched for the
        minimal set of constants. It is best to provide these for large system
        equations, as the search can be prohibitively long in duration.
    specifieds_symbols : iterable of Functions of time, optional
        If provided, the system's equations will not be searched for the
        minimal set of specifieds. It is best to provide these for large system
        equations, as the search can be prohibitively long in duration.

    """
    def __init__(self, eom_method, constants=None, specifieds=None,
                 ode_solver=None, initial_conditions=None, times=None,
                 outputs=None, noncontributing_forces=None,
                 constants_symbols=None, specifieds_symbols=None):

        self._last_generated_ode_user_kwargs = {}

        self._eom_method = eom_method
        # TODO : What if user adds symbols after constructing a System?
        if constants_symbols is None:
            self._constants_symbols = self._Kane_constant_symbols()
        else:
            self._constants_symbols = set(constants_symbols)
        if specifieds_symbols is None:
            self._specifieds_symbols = self._Kane_undefined_dynamicsymbols()
        else:
            self._specifieds_symbols = set(specifieds_symbols)

        self._extract_constraints()

        self._auxiliaries = []

        if outputs is None:
            outputs = dict()

        if self.constraints:
            self._con_syms = tuple(sm.Dummy('c' + str(i)) for i in
                                   range(self.num_constraints))
            outputs[self._con_syms] = self.constraints

        self.outputs = outputs  # calls _parse_outputs

        # NOTE: must be set before the state variables are intialized, so do
        # this first
        if noncontributing_forces is None:
            self._noncontributing_forces = []
        else:
            # calls parse_outputs again:
            self.noncontributing_forces = list(noncontributing_forces)

        if constants is None:
            self.constants = dict()
        else:
            self.constants = constants

        if specifieds is None:
            self.specifieds = dict()
        else:
            self.specifieds = specifieds

        if ode_solver is None:
            self.ode_solver = odeint
        else:
            self.ode_solver = ode_solver

        if initial_conditions is None:
            self.initial_conditions = dict()
        else:
            self.initial_conditions = initial_conditions

        if times is None:
            self.times = []  # gets converted to empty array([])
        else:
            self.times = times

        self._evaluate_ode_function = None
        self._needs_code_regeneration = True

    @property
    def coordinates(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized coordinates."""
        return self.eom_method.q[:]

    @property
    def num_coordinates(self):
        """Returns the number of coordinates."""
        return len(self.coordinates)

    @property
    def speeds(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized speeds."""
        return self.eom_method.u[:]

    @property
    def num_speeds(self):
        """Returns the number of speeds."""
        return len(self.speeds)

    @property
    def auxiliaries(self):
        """Returns a list of the symbols representing the system's auxiliary
        states which are the time integrals of any outputs that are linear
        functions of the time derivatices of the generalized speeds."""
        return self._auxiliaries

    @property
    def num_auxiliaries(self):
        """Returns the number of auxiliaries."""
        return len(self.auxiliaries)

    @property
    def states(self):
        """Returns a list of the symbolic functions of time representing the
        system's states, i.e. generalized coordinates plus the generalized
        speeds. These are in the same order as used in integration (as
        passed into evaluate_ode_function) and match the order of the mass
        matrix and forcing vector.

        """
        # requires settable attributes : outputs
        if self._linear_outputs_symbols:
            speeds = self.speeds + self.auxiliaries
            return self.coordinates + speeds
        else:
            return self.coordinates + self.speeds

    @property
    def num_states(self):
        """Returns the number of states."""
        return len(self.states)

    @property
    def eom_method(self):
        """This is a
        :external+sympy:py:class:`~sympy.physics.mechanics.kane.KanesMethod`.
        The method used to generate the equations of motion. Read-only."""
        return self._eom_method

    @property
    def constants(self):
        """A dict that provides the numerical values for the constants in the
        problem (all non-dynamics symbols). Keys are the symbols for the
        constants, and values are floats. Constants that are not specified in
        this dict are given a default value of 1.0.

        """
        return self._constants

    @constants.setter
    def constants(self, constants):
        self._check_constants(constants)
        self._constants = constants

    @property
    def num_constants(self):
        """Returns the number of constants."""
        return len(self.constants_symbols)

    @property
    def constants_symbols(self):
        """A set of the symbolic constants (not functions of time) in the
        system.
        """
        # requires settable attributes : constants
        return self._constants_symbols

    def _check_constants(self, constants):
        symbols = self.constants_symbols
        for k in constants.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a constant.".format(k))

    def _constants_padded_with_defaults(self):
        d = dict(zip(self.constants_symbols, repeat(1.0, self.num_constants)))
        d.update(self.constants)
        return d

    @property
    def _constants_array(self):
        p_dict = self._constants_padded_with_defaults()
        return np.array([p_dict[pi] for pi in self.constants_symbols])

    @property
    def specifieds(self):
        """A dict that provides numerical values for the specified quantities
        in the problem (all dynamicsymbols that are not defined by the
        equations of motion). There are two possible formats. (1) is more
        flexible, but (2) is more efficient (by a factor of 3).

        (1) Keys are the symbols for the specified quantities, or a tuple of
        symbols, and values are the floats, arrays of floats, or functions that
        generate the values. If a dictionary value is a function, it must have
        the same signature as ``f(x, t)``, the ode right-hand-side function
        (see the documentation for the :py:attr:`ode_solver` attribute). You
        needn't provide values for all specified symbols. Those for which you
        do not give a value will default to 0.0.

        (2) There are two keys: 'symbols' and 'values'. The value for 'symbols'
        is an iterable of *all* the specified quantities in the order that you
        have provided them in 'values'. Values is an ndarray, whose length is
        :py:attr:`num_specifieds`, or a function of x and t that returns an
        ndarray (also of length :py:attr:`num_specifieds`). NOTE: You must
        provide values for all specified symbols. In this case, we do *not*
        provide default values.

        NOTE: If you switch formats with the same instance of System, you
        *must* call :py:meth:`~pydy.system.System.generate_ode_function` before
        calling :py:meth:`~pydy.system.System.integrate` again.

        Examples
        --------
        Here are examples for (1). Keys can be individual symbols, or a tuple
        of symbols. Length of a value must match the length of the
        corresponding key. Values can be functions that return iterables::

            sys = System(km)
            sys.specifieds = {(a, b, c): np.ones(3), d: lambda x, t: -3 * x[0]}
            sys.specifieds = {(a, b, c): lambda x, t: np.ones(3)}

        Here are examples for (2)::

            sys.specifieds = {'symbols': (a, b, c, d),
                              'values': np.ones(4)}
            sys.specifieds = {'symbols': (a, b, c, d),
                              'values': lambda x, t: np.ones(4)}

        """
        return self._specifieds

    @specifieds.setter
    def specifieds(self, specifieds):
        self._check_specifieds(specifieds)
        self._specifieds = specifieds

    @property
    def num_specifieds(self):
        """Returns the number of specifieds."""
        return len(self.specifieds_symbols)

    @property
    def specifieds_symbols(self):
        """A set of the dynamicsymbols you must specify."""
        # TODO : Eventually use a method in the KanesMethod class.
        return self._specifieds_symbols

    def _assert_is_specified_symbol(self, symbol, all_symbols):
        if symbol not in all_symbols:
            raise ValueError("Symbol {} is not a 'specified' symbol.".format(
                symbol))

    def _assert_symbol_appears_multiple_times(self, symbol, symbols_so_far):
        if symbol in symbols_so_far:
            raise ValueError("Symbol {} appears more than once.".format(
                symbol))

    def _specifieds_are_in_format_2(self, specifieds):
        keys = specifieds.keys()
        if ('symbols' in keys and 'values' in keys):
            return True
        else:
            return False

    def _check_specifieds(self, specifieds):
        symbols = self.specifieds_symbols

        symbols_so_far = list()

        if self._specifieds_are_in_format_2(specifieds):

            # The symbols must be specifieds.
            for sym in specifieds['symbols']:
                self._assert_is_specified_symbol(sym, symbols)

            # Each specified symbol can appear only once.
            for sym in specifieds['symbols']:
                self._assert_symbol_appears_multiple_times(sym, symbols_so_far)
                symbols_so_far.append(sym)

            # Must have provided all specifieds.
            for sym in self.specifieds_symbols:
                if sym not in specifieds['symbols']:
                    raise ValueError(
                        "Specified symbol {} is not provided.".format(sym))

        else:

            for k, v in specifieds.items():

                # The symbols must be specifieds.
                if isinstance(k, tuple):
                    for ki in k:
                        self._assert_is_specified_symbol(ki, symbols)
                else:
                    self._assert_is_specified_symbol(k, symbols)

                # Each specified symbol can appear only once.
                if isinstance(k, tuple):
                    for ki in k:
                        self._assert_symbol_appears_multiple_times(
                            ki, symbols_so_far)
                        symbols_so_far.append(ki)
                else:
                    self._assert_symbol_appears_multiple_times(
                        k, symbols_so_far)
                    symbols_so_far.append(k)

    def _symbol_is_in_specifieds_dict(self, symbol, specifieds_dict):
        for k in specifieds_dict.keys():
            if symbol == k or (isinstance(k, tuple) and symbol in k):
                return True
        return False

    def _specifieds_padded_with_defaults(self):
        d = dict(zip(self.specifieds_symbols,
                     repeat(0.0, self.num_specifieds)))
        d.update(self.specifieds)
        return d

    @property
    def times(self):
        """A 1D ndarray of monotonic time values over which the equations of
        motion are numerically integrated. Can be set with an array-like for a
        shape(n,) array."""
        return self._times

    @times.setter
    def times(self, new_times):
        times = np.asarray(new_times)
        assert self._check_times(times)
        self._times = times

    def _check_times(self, times):

        # TODO : this check can probably be removed.
        if len(times.shape) == 0:
            raise TypeError("Times should be in an array_like format.")

        if not np.all(times >= 0):
            raise ValueError("Times supplied must have positive values.")

        if not np.all(np.diff(times) >= 0):
            raise ValueError("Times supplied should be in an ascending order.")

        return True

    @property
    def ode_solver(self):
        """A function that performs forward integration. It must have the same
        signature as :external+scipy:py:func:`scipy.integrate.odeint`, which
        is::

            x_history = ode_solver(f, x0, t, args=f_args)

        where ``f`` is a function ``f(x, t, *f_args)``, ``x0`` are the initial
        conditions, ``x_history`` is the state time history, ``x`` is the
        state, ``t`` is the time, and ``args`` is a keyword argument takes
        arguments that are then passed to ``f``. The default solver is
        :external+scipy:py:func:`scipy.integrate.odeint`.

        Examples
        ========

        SciPy introduced a unified :py:func:`scipy.integrate.solve_ivp` API
        which can be used with PyDy. ``solve_ivp`` requires a function that has
        swapped first arguments and it returns a solution object where the
        trajectory is the transpose of what ``odeint`` outputs. You can make a
        custom ODE solver function to use ``solve_ivp`` like so:

        >>> from pydy.models import multi_mass_spring_damper
        >>> sys = multi_mass_spring_damper()
        >>> sys.initial_conditions[sys.coordinates[0]] = 1.0
        >>> sys.times = [1.0, 2.0, 3.0]
        >>> from scipy.integrate import solve_ivp
        >>> def custom_ode_solver(f, x0, ts, args=(), **kwargs):
        ...     return solve_ivp(lambda t, x: f(x, t, *args), ts[[0, -1]], x0,
        ...                      t_eval=ts, **kwargs).y.T
        >>> sys.ode_solver = custom_ode_solver

        This then allows one to easiliy change methods and settings following
        SciPy's API:

        >>> sys.integrate(method='LSODA', rtol=1e-10)
        array([[ 1.00000000e+00, -5.67952532e-17],
               [ 6.59700039e-01, -5.33506568e-01],
               [ 1.50574778e-01, -4.19279930e-01]])
        >>> sys.integrate(method='RK23', rtol=1e-12)
        array([[ 1.        ,  0.        ],
               [ 0.65970115, -0.53350624],
               [ 0.15057689, -0.41928088]])

        """
        return self._ode_solver

    @ode_solver.setter
    def ode_solver(self, ode_solver):
        if not hasattr(ode_solver, '__call__'):
            msg = "``ode_solver`` ({}) is not a function."
            raise ValueError(msg.format(ode_solver))
        self._ode_solver = ode_solver

        # NOTE : This ensures that force_c_contiguous will be set to True if a
        # cutom integrator is added and you last generated a cython based ode
        # function.
        if 'generator' in self._last_generated_ode_user_kwargs:
            if (self._last_generated_ode_user_kwargs['generator'] == 'cython'
                    and ode_solver is not odeint):
                self._needs_code_regeneration = True

    @property
    def initial_conditions(self):
        """Initial conditions for all states (coordinates and speeds). Keys
        are the symbols for the coordinates and speeds, and values are
        floats. Coordinates or speeds that are not specified in this dict
        are given a default value of 0.0.

        """
        return self._initial_conditions

    @initial_conditions.setter
    def initial_conditions(self, initial_conditions):
        self._check_initial_conditions(initial_conditions)
        self._initial_conditions = initial_conditions

    def _check_initial_conditions(self, initial_conditions):
        symbols = self.states
        for k in initial_conditions.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a state.".format(k))

    def _initial_conditions_padded_with_defaults(self):
        d = dict(zip(self.states, repeat(0.0, self.num_states)))
        d.update(self.initial_conditions)
        return d

    @property
    def _initial_conditions_array(self):
        x0_dict = self._initial_conditions_padded_with_defaults()
        return np.array([x0_dict[xi] for xi in self.states])

    @property
    def outputs(self):
        """Dictionary of functions of time or utple of functions of time mapped
        to SymPy expressions or iterables of expressions that represent extra
        functions of the state that should be evaluted alongside the ordinary
        differential equations. Acceptable key pairs for this dictionary take
        the following three forms:

        A single function of time mapped to a function of the state::

            outputs[p(t)] = k*x(t)

        A tuple of functions of time mapped to functions of the state::

            outputs[(f1(t), f2(t))] = (k*x(t), c*v(t))

        A tuple of functions of time mapped to a system of linear equations in
        the functions and the time derivatives of the states::

            outputs[(m1(t), m2(t))] = (m1(t) - 4*m2(t) + k*v(t).diff(t) + 2,
                                       m1(t) + 3*m2(t) - omega(t).diff(t))

        If equations of the last form are provided, this linear system will be
        numerically solved alongside the ordinary differential equations.

        Notes
        =====

        If your system has configuration or motion constraints, these will
        automatically be added to the outputs dictionary. If your system has
        noncontributing forces exposed and you provide names for those forces,
        these will automatically be added to the outputs dictionary.

        """
        return self._outputs

    @outputs.setter
    def outputs(self, outputs):
        self._outputs = outputs
        self._parse_outputs()
        # NOTE : It the output equations are updated they may have new symbols.
        exprs = []
        if self._simple_outputs_symbols:
            exprs += self._simple_outputs_matrix[:]
        if self._linear_outputs_symbols:
            exprs += self._linear_outputs_mass_matrix_rows[:]
            exprs += self._linear_outputs_forcing_rows[:]
        for s in self._Kane_constant_symbols(exprs=exprs):
            self._constants_symbols.add(s)
        for s in self._Kane_undefined_dynamicsymbols(exprs=exprs):
            self._specifieds_symbols.add(s)
        self._needs_code_regeneration = True

    def _parse_outputs(self):
        # Divide the equations into three types:
        #
        # 1. Equations that are simply a function of the state:
        # [y1] = [f1(t, x, r, p)]
        # [y2]   [f2(t, x, r, p)]
        # [y3]   [f3(t, x, r, p)]
        #
        # 2. Equations that are linear in the state derivatives and the new
        # output variables, for example noncontributing forces. The essential
        # equations of motion will be augmented with these equations and the
        # outputs y will be solved for along in the inversion of the new
        # augmented mass matrix.
        # [Md  0] [u'] = [Fd]
        # [Mu My] [y ]   [Fy]
        #
        # TODO : How could we also support Lagrange multipliers?:
        # [Md CT] [u'] = [Fd]
        # [C  Ml] [l ]   [Fl]
        #
        # TODO : support nonlinear functions of x' later.
        # 3. Equations that are nonlinear functions of the state and its state
        # derivatives.
        # [y1] = [f1(t, x', x, r, p)]
        # [y2]   [f2(t, x', x, r, p)]
        # [y3]   [f3(t, x', x, r, p)]
        #
        # The end goal is to have something like:
        #
        # def rhs(t, x, r, p):
        #     M, F, Y1 = eval_eqs(t, x, r, p)
        #     sol = solve(M, F)
        #     xdot = sol[:len(x)]
        #     Y2 = sol[len(x):]
        #     Y3 = eval_eqs2(t, x, xdot, r, p)
        #     Y = put_in_order((Y1, Y2, Y3))
        #     return xdot, Y
        #
        # We should retain the order of the outputs in the provided dictionary
        # for Y.
        #
        # The dictionary is iterated and for all ouputs:
        # Y = [y1, y2, y3, y4, y5, y6]
        # it is separated into simple outputs:
        # Y1 = [y1, y2, y3, y6]
        # and linear outputs:
        # Y2 = [y4, y5]
        # so to reconstruct Y in the order of the dictionary we need to store
        # the indices in Y so we can do:
        # Y[simple_idxs] = Y1
        # Y[linear_idxs] = Y2

        output_names_in_order = []

        funcs_of_x = []
        simple_outputs_names = []

        funcs_of_xdot = []
        linear_eq_names = []

        # TODO : KanesMethod should store the linear components of the
        # auxiliary equations. It should also have a method/attribute to return
        # the augmented mass matrix and forcing vector.

        for var, expr in self.outputs.items():
            if isinstance(var, tuple):
                for v, e in zip(var, expr):
                    output_names_in_order.append(v)
                    if e.has(sm.Derivative):
                        funcs_of_xdot.append(e)
                        linear_eq_names.append(v)
                    else:
                        funcs_of_x.append(e)
                        simple_outputs_names.append(v)
            else:
                output_names_in_order.append(var)
                if expr.has(sm.Derivative):
                    funcs_of_xdot.append(expr)
                    linear_eq_names.append(var)
                else:
                    funcs_of_x.append(expr)
                    simple_outputs_names.append(var)

        if len(set(output_names_in_order)) < len(output_names_in_order):
            raise ValueError('All outputs must have unique names.')

        self._num_simple_outputs = len(simple_outputs_names)
        self._simple_outputs_symbols = simple_outputs_names
        if funcs_of_x:
            self._simple_outputs_matrix = sm.Matrix(funcs_of_x)
        else:
            self._simple_outputs_matrix = funcs_of_x

        if funcs_of_xdot:
            funcs_of_xdot = sm.Matrix(funcs_of_xdot)
            xd = [ui.diff() for ui in self.speeds] + linear_eq_names
            mass_matrix_rows, forcing_rows = sm.linear_eq_to_matrix(
                funcs_of_xdot, xd)
        else:
            mass_matrix_rows = sm.Matrix([])
            forcing_rows = sm.Matrix([])

        self._auxiliaries = [sm.Symbol('∫ ' + s.name + ' dt') for s in
                             linear_eq_names]
        self._num_linear_outputs = len(linear_eq_names)
        self._linear_outputs_symbols = linear_eq_names
        self._linear_outputs_mass_matrix_rows = mass_matrix_rows
        self._linear_outputs_forcing_rows = forcing_rows

        if self.constraints:
            self._constraint_idxs = [
                self._simple_outputs_symbols.index(ci) for ci in
                self._con_syms]

        self.outputs_symbols = output_names_in_order
        self._num_outputs = len(output_names_in_order)

        self._simple_idxs = [output_names_in_order.index(si)
                             for si in simple_outputs_names]
        self._linear_idxs = [output_names_in_order.index(si)
                             for si in linear_eq_names]

    @property
    def num_outputs(self):
        """Returns the number of outputs."""
        return self._num_outputs

    @property
    def noncontributing_forces(self):
        """List of symbolic functions of time representing the noncontributing
        forces (force & torque measure numbers) associated with auxiliary
        speeds."""
        return self._noncontributing_forces

    @noncontributing_forces.setter
    def noncontributing_forces(self, noncontributing_forces):
        if not hasattr(self.eom_method, 'auxiliary_eqs'):
            msg = ('The KanesMethod object has no auxiliary equations and '
                   'thus noncontributing forces cannot be provided.')
            raise RuntimeError(msg)
        if len(noncontributing_forces) != len(self.eom_method._uaux):
            msg = ('You must provide symbols for {} noncontributing forces '
                   'that are present in the auxiliary equations.')
            raise ValueError(msg.format(len(self.eom_method._uaux)))
        # TODO : Check that the noncontributing force symbols are present in
        # the auxiliary equations and that they are not a coordinate, speed, or
        # specified.
        # TODO : What is this check?
        self._noncontributing_forces = list(noncontributing_forces)
        if tuple(noncontributing_forces) in self.outputs:
            raise ValueError('Constraint loads already present in outputs.')

        non_syms = tuple(noncontributing_forces)
        self.outputs[non_syms] = self.eom_method.auxiliary_eqs
        self._parse_outputs()
        self._needs_code_regeneration = True

    @property
    def evaluate_ode_function(self):
        """A function generated by
        :py:func:`~pydy.codegen.ode_function_generators.generate_ode_function`
        that computes the state derivatives::

            xd = evaluate_ode_function(x, t, *args)

        This function is used by the :py:attr:`~pydy.system.System.ode_solver`.

        To see the autogenerated docstring and expected arguments call
        :py:func:`help`::

            help(system.evaluate_ode_function)

        """
        # TODO : It would be more useful if the generated docstring was shown
        # when system.evaluate_ode_function? in IPython is called.
        # Interestingly help(system.evaluate_ode_function) does show the
        # underlying function's docstring.
        return self._evaluate_ode_function

    def _args_for_gen_ode_func(self):
        """Returns a tuple of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generate_ode_function``.

        """
        if self._linear_outputs_symbols:
            Fd = self.eom_method.forcing
            Fa = self._linear_outputs_forcing_rows
            # [Md  0] [u'] = [Fd]
            # [Mu Mj] [j']   [Fa]
            forcing = Fd.col_join(Fa)
            speeds = self.speeds + self.auxiliaries
        else:
            forcing = self.eom_method.forcing
            speeds = self.speeds

        args = (forcing,
                self.coordinates,
                speeds,
                list(sm.ordered(self.constants_symbols)))

        return args

    def _kwargs_for_gen_ode_func(self):
        """Returns a dictionary of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generage_ode_function``.

        """
        if self._specifieds_are_in_format_2(self.specifieds):
            specifieds = self.specifieds['symbols']
        else:
            specifieds = self.specifieds_symbols

        # generate_ode_func does not accept an empty tuple for the
        # specifieds, so set it to None
        if not specifieds:
            specifieds = None

        kin_diff_dict = self.eom_method.kindiffdict()
        kin_diff_rhs = sm.Matrix([kin_diff_dict[q.diff()] for q in
                                  self.coordinates])

        if self._linear_outputs_symbols:
            Md = self.eom_method.mass_matrix
            MuMj = self._linear_outputs_mass_matrix_rows
            Mz = sm.zeros(Md.shape[0], MuMj.shape[0])
            # [Md  0] [u'] = [Fd]
            # [Mu Mj] [j']   [Fa]
            mass_matrix = Md.row_join(Mz).col_join(MuMj)
        else:
            mass_matrix = self.eom_method.mass_matrix

        kwargs = {
            'mass_matrix': mass_matrix,
            'coordinate_derivatives': kin_diff_rhs,
            'specifieds': specifieds,
        }

        if self._simple_outputs_symbols:
            kwargs['outputs'] = self._simple_outputs_matrix

        return kwargs

    def _extract_constraints(self):
        """Extracts the configuration, nonholonomic, and motion constraints
        from the eom_method and stores them in attributes."""

        if self.eom_method._f_h:
            self.holonomic_constraints = self.eom_method._f_h
        else:
            self.holonomic_constraints = sm.Matrix([])

        self._num_holonomic_constraints = len(self.holonomic_constraints)

        if self.eom_method._k_nh:
            # TODO : KanesMethod and _Method should store the original
            # constraints passed by the user. Fix in sympy.physics.mechanics!

            # rebuild the velocity constraints from KanesMethod, note that
            # these can include time differentiated holonomic constraints
            velocity_constraints = (
                self.eom_method._k_nh*self.eom_method.u +
                self.eom_method._f_nh)

            self.velocity_constraints = velocity_constraints

            if self.holonomic_constraints:
                h_idxs, nh_idxs = _sort_velocity_constraints(
                    velocity_constraints,
                    self.coordinates[:],
                    list(self.eom_method.kindiffdict().values()))

                if len(h_idxs) != self.num_holonomic_constraints:
                    raise ValueError('There should be the same number of time'
                                    ' differentiated holonomic constraints '
                                    ' present in the velocity constraints to '
                                    ' that of the holonomic constraints.')
                if len(nh_idxs) == 0:
                    self.nonholonomic_constraints = sm.Matrix([])
                else:
                    self.nonholonomic_constraints = sm.Matrix(
                        velocity_constraints)[nh_idxs, 0]
            else:
                self.nonholonomic_constraints = velocity_constraints
        else:
            self.velocity_constraints = sm.Matrix([])
            self.nonholonomic_constraints = sm.Matrix([])

        self._num_nonholonomic_constraints = len(self.nonholonomic_constraints)
        self._num_velocity_constraints = len(self.velocity_constraints)

    @property
    def num_holonomic_constraints(self):
        """Number of configuration constraints."""
        return self._num_holonomic_constraints

    @property
    def num_nonholonomic_constraints(self):
        """Number of nonholonomic constraints."""
        return self._num_nonholonomic_constraints

    @property
    def num_velocity_constraints(self):
        """Number of motion constraints."""
        return self._num_velocity_constraints

    @property
    def constraints(self):
        """A column matrix of configuration and nonholonomic constraints
        expressions, ordered as stored in
        :external+sympy:py:class:`~sympy.physics.mechanics.kane.KanesMethod`."""
        constraints = sm.Matrix([])

        if self.holonomic_constraints or self.nonholonomic_constraints:
            if self.holonomic_constraints:
                constraints = self.holonomic_constraints

            if constraints and self.nonholonomic_constraints:
                constraints = constraints.col_join(
                    self.nonholonomic_constraints)
            else:
                constraints = self.nonholonomic_constraints

        return constraints

    @property
    def num_constraints(self):
        """Total number of configuration and nonholonomic constaints."""
        return (self.num_holonomic_constraints +
                self.num_nonholonomic_constraints)

    def set_dependent_initial_conditions(self, dep_vars=None, use_jac=False,
                                         **root_kwargs):
        """Sets the initial conditions of the dependent coordinates and
        dependent speeds using the holonomic and nonholonomic constraints,
        respectively.

        Parameters
        ==========
        dep_vars : iterable of Function()(t), optional
            Dependent coordinates and speeds to solve for. The number of
            coordinates should be equal to the number of holonomic constraints.
            The number of speeds should be equal to the number of nonholonic
            constraints. If None, the dependent coordinates and speeds are
            those used in KanesMethod instantiation.
        use_jac : boolean, optional
            If true the Jacobian of the constraint equations will be used to
            solve the constraint equations for the dependent states.
        root_kwargs
            Extra keyword arguments that are passed to
            :external+scipy:py:func:`scipy.optimize.root`.

        """
        # TODO : The nonholonomic constraints can be solved analytically, root
        # is only required for the coordinates.

        if self.num_constraints == 0:
            msg = ('This system does not have constraints, set all initial '
                   'conditions yourself.')
            raise ValueError(msg)
        else:
            num_holo = self.num_holonomic_constraints
            num_nonh = self.num_nonholonomic_constraints

        # TODO : These variables should be publicly accessible on KanesMethod.
        # NOTE : If there are holonomic constraints, then you have to solve for
        # both the dependent coordinate and dependent speed associated with
        # this constraint, i.e. the time differentiated holoomic constraints is
        # treated as a nonholonomic constraint.
        if dep_vars is None:
            dep_vars = self.eom_method._qdep[:] + self.eom_method._udep[:]

        # TODO : Would be nice to check if the dependent variables are present
        # in the constraints and that the right number of coordinates and
        # speeds are each supplied.
        if len(dep_vars) != (2*num_holo + num_nonh):
            msg = (f'You must supply {num_holo} dependent coordinates and '
                   f'{num_nonh} dependent speeds.')
            raise ValueError(msg)

        x = self._initial_conditions_array
        p = self._constants_array

        x0_dict = self._initial_conditions_padded_with_defaults()
        dep_guess = [x0_dict[xi] for xi in dep_vars]
        dep_idxs = [self.states.index(xi) for xi in dep_vars]

        # NOTE : M + (M + m) equations.
        # TODO : Create a ._evaluate_all_constraints() method for this.
        constraints = self.holonomic_constraints.col_join(
            self.velocity_constraints)

        if use_jac:
            jac = constraints.jacobian(dep_vars)
            eval_jac = sm.lambdify((self.states, self.constants_symbols),
                                   (constraints, jac), cse=True)

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                con_vals, jac_vals = eval_jac(x, p)
                return con_vals.squeeze(), jac_vals

            fprime = True
        else:
            # TODO : Create a ._evaluate_all_constraints() method for this.
            eval_con = sm.lambdify((self.states, self.constants_symbols),
                                   constraints, cse=True)

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                return eval_con(x, p).squeeze()

            fprime = False

        # not settable by the user, so overwrite:
        root_kwargs['args'] = (p, )
        root_kwargs['jac'] = fprime

        sol = root(eval_f, dep_guess, **root_kwargs)
        if not sol.success:
            msg = ('Failed to find a solution that meets tolerance. Maybe a '
                   'better guess will help or you may have to manually solve '
                   'for the dependent coordinates. SciPy root() failure '
                   'message: ' + sol.message)
            warnings.warn(msg, PyDyUserWarning, stacklevel=2)

        dep_vals = sol.x

        for si, vi in zip(dep_vars, dep_vals):
            self.initial_conditions[si] = vi

    def generate_ode_function(self, **kwargs):
        """Returns a function generated from
        :py:func:`~pydy.codegen.ode_function_generators.generate_ode_function`
        with the appropriate arguments and also sets the
        ``evaluate_ode_function`` attribute to the resulting function.

        Parameters
        ----------
        kwargs
            All other kwargs are passed onto
            :py:func:`pydy.codegen.ode_function_generators.generate_ode_function`.
            Don't specify the ``specifieds`` keyword argument though; the
            ``System`` class takes care of those.

        Returns
        -------
        evaluate_ode_function : function
            A function which evaluates the derivaties of the states.

        Notes
        -----

        If the Cython generator is selected and you have a custom
        ``ode_solver`` set, keyword argument ``force_c_contiguous`` will be
        automatically set to ``True``. You can disable this by setting it to
        ``False`` but you must ensure ensure that ode solver only passes C
        contiguous arrays to the generated ode function. Forcing C contiguous
        arrays introduces a small performance penalty due to the necessity of
        copying arrays.

        """
        self._parse_outputs()  # call before generating the args/kwargs below

        self._last_generated_ode_user_kwargs = kwargs.copy()

        if 'specified' in kwargs:
            kwargs.pop('specified')
            print("User supplied 'specified' kwarg was disregarded.")

        if 'specifieds' in kwargs:
            kwargs.pop('specifieds')
            print("User supplied 'specifieds' kwarg was disregarded.")

        kwargs.update(self._kwargs_for_gen_ode_func())

        if 'force_c_contiguous' not in kwargs:
            if 'generator' in kwargs:
                if (kwargs['generator'] == 'cython' and self.ode_solver is
                        odeint):
                    kwargs['force_c_contiguous'] = False
                # NOTE : This ensures that the arrays are forced to be C
                # contiguous if a user sets an integrator other than odeint,
                # but only if they are using the Cython generator. There is a
                # performance penalty but it avoids errors and they can
                # manually override this to False if they know their ode_solver
                # choice delivers C contiguous arrays.
                elif kwargs['generator'] == 'cython':
                    kwargs['force_c_contiguous'] = True

        self._evaluate_ode_function = generate_ode_function(
            *self._args_for_gen_ode_func(),
            **kwargs)

        return self.evaluate_ode_function

    def _prep_for_evaluate(self):
        # Users might have changed these properties by directly accessing the
        # dict, without using the setter. Before we integrate, make sure they
        # did not muck up these dicts.
        self._check_constants(self.constants)
        self._check_specifieds(self.specifieds)
        self._check_initial_conditions(self.initial_conditions)
        assert self._check_times(self.times)

        if self.evaluate_ode_function is None or self._needs_code_regeneration:
            self.generate_ode_function(**self._last_generated_ode_user_kwargs)
            self._needs_code_regeneration = False

        if self._specifieds_are_in_format_2(self.specifieds):
            specified_value = self.specifieds['values']
        else:
            specified_value = self._specifieds_padded_with_defaults()

        # If there are no specifieds then specified_value will be an empty
        # dict.
        if isinstance(specified_value, dict) and not specified_value:
            args = (self._constants_padded_with_defaults(),)
        else:
            args = (specified_value, self._constants_padded_with_defaults())

        return self._initial_conditions_array, args

    def _prep_x_t_overrides(self, x, t):
        x_default, args = self._prep_for_evaluate()

        if x is None:
            x = x_default
        x = np.asarray(x)

        # if times has not been set, just use t=0.0
        if t is None:
            if len(x.shape) == 1:  # x at t
                if self.times.size == 0:  # array([])
                    t = 0.0
                else:
                    t = self.times[0]
            else:
                if self.times.size == 0:  # array([])
                    t = np.zeros(x.shape[0])
                else:
                    t = self.times

        return x, t, args

    def evaluate_ode(self, x=None, t=None):
        """Returns the right hand side of the differential equations. The
        default is to evaluate at the set initial_conditions at the first time
        value or with t=0 if :py:attr:`times` is not set. Pass in optional
        arguments to override using the initial state and time.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State values at time t.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        xd : ndarray, shape(n,) or shape(m, n)
           Time derivative of the states at time t.

        Notes
        =====

        This method is present for convenience, it is not designed to be used
        where performance matters, use
        :py:attr:`~pydy.system.System.evaluate_ode_function` directly when
        performance is needed.

        To see the order of the state values use::

            system = System(...)
            system.states

        or::

            rhs = system.generate_ode_function()
            help(rhs)

        """
        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            res = self.evaluate_ode_function(x, t, *args)
            if self._simple_outputs_symbols:
                return res[0]
            else:
                return res
        # NOTE : I tried to make use of numpy.vectorize but it is not possible
        # due to args not being necessarily being comprised of arrays.
        elif len(x.shape) == 2:
            if isinstance(t, float):
                raise ValueError('t must be array with the same length as x.')
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            xd = np.zeros_like(x)
            for i, (ti, xi) in enumerate(zip(t, x)):
                res = self.evaluate_ode_function(xi, ti, *args)
                if self._simple_outputs_symbols:
                    xd[i, :] = res[0]
                else:
                    xd[i, :] = res
            return xd

    def evaluate_outputs(self, x=None, t=None):
        """Returns an array of the evaluated outputs. The default is to
        evaluate at the initial conditions at the first time value. Pass in
        optional arguments to override the state or time.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State values at time t.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        y : ndarray, shape(o,) or shape(m, o)
           o output values at time t.

        Notes
        =====

        This method is present for convenience, it is not designed to be used
        where performance matters, use
        :py:attr:`~pydy.system.System.evaluate_ode_function` directly when
        performance is needed.

        To see the order of the state values use::

            system = System(...)
            system.states

        or::

            rhs = system.generate_ode_function()
            help(rhs)

        """
        if not self.outputs:
            raise ValueError('This system has no outputs.')

        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            if self._linear_outputs_symbols and self._simple_outputs_symbols:
                y = np.zeros(self.num_outputs)
                xdot, y1 = self.evaluate_ode_function(x, t, *args)
                y[self._simple_idxs] = y1
                y[self._linear_idxs] = xdot[-len(self.auxiliaries):]
                return y
            elif (self._linear_outputs_symbols and not
                  self._simple_outputs_symbols):
                xdot = self.evaluate_ode_function(x, t, *args)
                return xdot[-len(self.auxiliaries):]
            else:
                return self.evaluate_ode_function(x, t, *args)[1]
        # NOTE : I tried to make use of numpy.vectorize but it is not possible
        # due to args not being necessarily being comprised of arrays.
        elif len(x.shape) == 2:
            if isinstance(t, float):
                raise ValueError('t must be array with the same length as x.')
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            y = np.zeros((len(t), self.num_outputs))
            for i, (ti, xi) in enumerate(zip(t, x)):
                if (self._linear_outputs_symbols and
                        self._simple_outputs_symbols):
                    xdot, y1 = self.evaluate_ode_function(xi, ti, *args)
                    y[i, self._simple_idxs] = y1
                    y[i, self._linear_idxs] = xdot[-len(self.auxiliaries):]
                elif (self._linear_outputs_symbols and not
                      self._simple_outputs_symbols):
                    xdot = self.evaluate_ode_function(xi, ti, *args)
                    y[i, :] = xdot[-len(self.auxiliaries):]
                else:
                    y[i, :] = self.evaluate_ode_function(xi, ti, *args)[1]
            return y

    def evaluate_constraints(self, x=None, t=None):
        """Returns the values of the configuration and motion constraints at
        the initial condition or, alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            system = System(...)
            system.states

        or::

            rhs = system.generate_ode_function()
            help(rhs)

        """
        # TODO : convert this to calling evaluate_outputs() and then selecting
        # the constraint values.
        if self.num_constraints == 0:
            raise ValueError('This system has no constraints.')

        x, t, args = self._prep_x_t_overrides(x, t)

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            y = self.evaluate_ode_function(x, t, *args)[1]
            return y[self._constraint_idxs]
        elif len(x.shape) == 2:
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            con = np.zeros((x.shape[0], self.num_constraints))
            for i, (ti, xi) in enumerate(zip(t, x)):
                y = self.evaluate_ode_function(xi, ti, *args)[1]
                con[i, :] = y[self._constraint_idxs]
            return con

    def evaluate_holonomic_constraints(self, x=None, t=None):
        """Returns the values of the configuration at the initial condition or,
        alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            system = System(...)
            system.states

        or::

            rhs = system.generate_ode_function()
            help(rhs)

        """
        if self.num_holonomic_constraints == 0:
            raise ValueError('This system has no configuration constraints.')
        con = self.evaluate_constraints(x=x, t=t)
        if len(con.shape) == 1:
            return con[:self.num_holonomic_constraints]
        else:
            return con[:, :self.num_holonomic_constraints]

    def evaluate_velocity_constraints(self, x=None, t=None):
        """Returns the values of the velocity constraints at the initial
        condition or, alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            system = System(...)
            system.states

        or::

            rhs = system.generate_ode_function()
            help(rhs)

        """
        if self.num_velocity_constraints == 0:
            raise ValueError('This system has no motion constraints.')
        con = self.evaluate_constraints(x=x, t=t)
        if len(con.shape) == 1:
            return con[self.num_holonomic_constraints:]
        else:
            return con[:, self.num_holonomic_constraints:]

    def integrate(self, **solver_kwargs):
        """Integrates the equations
        :py:attr:`~pydy.system.System.evaluate_ode_function` using
        :py:attr:`~pydy.system.System.ode_solver`.

        It is necessary to have first generated an ode function. If you have
        not done so, we do so automatically by invoking
        :py:meth:`~pydy.system.System.generate_ode_function`. However, if you
        want to customize how this function is generated (e.g., change the
        generator to cython), you can call
        :py:meth:`~pydy.system.System.generate_ode_function` on your own
        (before calling :py:meth:`~pydy.system.System.integrate`).

        Parameters
        ----------
        **solver_kwargs
            Optional arguments that are passed on to the
            :py:attr:`~pydy.system.System.ode_solver`.

        Returns
        -------
        x_history : ndarray, shape(num_integrator_time_steps, num_states)
            The trajectory of states (coordinates and speeds) through the
            requested time interval. num_integrator_time_steps is either
            len(times) if len(times) > 2, or is determined by the
            :py:attr:`~pydy.system.System.ode_solver`.

        """
        if len(self.times) < 2:
            raise ValueError('The times vector must be at least length 2.')

        x0, args = self._prep_for_evaluate()

        # NOTE : User cannot pass in args, System handles that.
        solver_kwargs.pop('args', None)

        # NOTE : skip the outputs if present
        def func(*args, **kwargs):
            return self.evaluate_ode_function(*args, **kwargs)[0]

        x_history = self.ode_solver(
            (func if self._simple_outputs_symbols else
             self.evaluate_ode_function),
            x0,
            self.times,
            args=args, **solver_kwargs)

        return x_history

    def _Kane_inlist_insyms(self):
        """TODO temporary."""

        uaux = self.eom_method._uaux[:]
        uauxdot = [sm.diff(i, dynamicsymbols._t) for i in uaux]

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.

        # TODO : KanesMethod should provide public attributes for qdot,
        # udot, uaux, and uauxdot.
        # NOTE : There have been breaking changes in SymPy, e.g.
        # https://github.com/pydy/pydy/issues/395, so all of these need to be
        # converted to lists before summing.
        insyms = set(list(self.eom_method.q[:]) +
                     list(self.eom_method._qdot[:]) +
                     list(self.eom_method.u[:]) +
                     list(self.eom_method._udot[:]) +
                     list(uaux) + list(uauxdot))

        inlist = (self.eom_method.forcing_full[:] +
                  self.eom_method.mass_matrix_full[:])

        return inlist, insyms

    def _Kane_undefined_dynamicsymbols(self, exprs=None):
        """Similar to ``_find_dynamicsymbols()``, except that it checks all
        syms used in the system. Code is copied from ``linearize()``.

        TODO temporarily here until KanesMethod and Lagranges method have an
        interface for obtaining these quantities.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        if exprs:
            from_eoms = exprs
        functions_of_time = set()
        for expr in from_eoms:
            functions_of_time = functions_of_time.union(
                find_dynamicsymbols(expr))
        return functions_of_time.difference(from_sym_lists)

    def _Kane_constant_symbols(self, exprs=None):
        """Similar to ``_find_othersymbols()``, except it checks all syms used
        in the system.

        Remove the time symbol.

        TODO temporary.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        if exprs:
            from_eoms = exprs
        unique_symbols = set()
        for expr in from_eoms:
            unique_symbols = unique_symbols.union(expr.free_symbols)
        constants = unique_symbols
        constants.remove(dynamicsymbols._t)
        return constants
