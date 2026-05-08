"""The System class manages the simulation (integration) of a system whose
equations are given by KanesMethod.

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

    km = KanesMethod(...)
    km.kanes_equations(force_list, body_list)
    times = np.linspace(0, 5, 100)
    sys = System(km, times=times)
    sys.integrate()

In this case, we use defaults for the numerical values of the constants,
specified quantities, initial conditions, etc. You probably won't like
these defaults. You can also specify such values via constructor keyword
arguments or via the attributes::

    sys = System(km,
                 initial_conditions={dynamicsymbol('q1'): 0.5},
                 times=times)
    sys.constants = {symbol('m'): 5.0}
    sys.integrate()

To double-check the constants, specifieds, states and times in your problem,
look at these properties::

    sys.constants_symbols
    sys.specifieds_symbols
    sys.states
    sys.times

In this case, the System generates the numerical ode function for you
behind the scenes. If you want to customize how this function is generated,
you must call ``generate_ode_function`` on your own::

    sys = System(KM)
    sys.generate_ode_function(generator='cython')
    sys.integrate()

"""
import warnings
from itertools import repeat

import numpy as np
import sympy as sm
from sympy.physics.mechanics import dynamicsymbols, find_dynamicsymbols
from scipy.integrate import odeint
from scipy.optimize import root

from .codegen.ode_function_generators import generate_ode_function
from .utils import PyDyFutureWarning, PyDyUserWarning

SYMPY_VERSION = sm.__version__


warnings.simplefilter('once', PyDyFutureWarning)


class System(object):
    """See the class's attributes for a description of the arguments to
    this constructor.

    The parameters to this constructor are all attributes of the System.
    Actually, they are properties. With the exception of ``eom_method``,
    these attributes can be modified directly at any future point.

    Parameters
    ----------
    eom_method : sympy.physics.mechanics.KanesMethod
        You must have called ``KanesMethod.kanes_equations()`` *before*
        constructing this ``System``.
    constants : dict, optional (default: all 1.0)
        This dictionary maps SymPy Symbol objects to floats.
    specifieds : dict, optional (default: all 0.0)
        This dictionary maps SymPy Functions of time objects, or tuples of
        them, to floats, NumPy arrays, or functions of the state and time.
    ode_solver : function, optional (default: scipy.integrate.odeint)
        This function computes the derivatives of the states.
    initial_conditions : dict, optional (default: all zero)
        This dictionary maps SymPy Functions of time objects to floats.
    times : array_like, shape(n,), optional
        An array_like object, which contains time values over which
        equations are integrated. It has to be supplied before
        System.integrate() can be called.

    """
    def __init__(self, eom_method, constants=None, specifieds=None,
                 ode_solver=None, initial_conditions=None, times=None):

        self._eom_method = eom_method

        # TODO : What if user adds symbols after constructing a System?
        # TODO : For large equations of motion, these two methods can take
        # unncessary time. One option is to let the user optionally pass in
        # these lists of symbols or derive them from the constants and
        # specifieds dictionaries.
        self._constants_symbols = self._Kane_constant_symbols()
        self._specifieds_symbols = self._Kane_undefined_dynamicsymbols()

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
            self.times = []
        else:
            self.times = times

        self._evaluate_ode_function = None

    @property
    def coordinates(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized coordinates."""
        return self.eom_method.q[:]

    @property
    def speeds(self):
        """Returns a list of the symbolic functions of time representing the
        system's generalized speeds."""
        return self.eom_method.u[:]

    @property
    def states(self):
        """Returns a list of the symbolic functions of time representing the
        system's states, i.e. generalized coordinates plus the generalized
        speeds. These are in the same order as used in integration (as
        passed into evaluate_ode_function) and match the order of the mass
        matrix and forcing vector.

        """
        return self.coordinates + self.speeds

    @property
    def eom_method(self):
        """This is a sympy.physics.mechanics.KanesMethod. The method used to
        generate the equations of motion. Read-only.

        """
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
    def constants_symbols(self):
        """A set of the symbolic constants (not functions of time) in the
        system.
        """
        return self._constants_symbols

    def _check_constants(self, constants):
        symbols = self.constants_symbols
        for k in constants.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a constant.".format(k))

    def _constants_padded_with_defaults(self):
        d = dict(zip(self.constants_symbols,
                     repeat(1.0, len(self.constants_symbols))))
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
        (see the documentation for the ``ode_solver`` attribute). You needn't
        provide values for all specified symbols. Those for which you do not
        give a value will default to 0.0.

        (2) There are two keys: 'symbols' and 'values'. The value for 'symbols'
        is an iterable of *all* the specified quantities in the order that you
        have provided them in 'values'. Values is an ndarray, whose length is
        `len(sys.specifieds_symbols)`, or a function of x and t that returns an
        ndarray (also of length `len(sys.specifieds_symbols)`). NOTE: You must
        provide values for all specified symbols. In this case, we do *not*
        provide default values.

        NOTE: If you switch formats with the same instance of System, you
        *must* call `generate_ode_function()` before calling `integrate()`
        again.

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
                     repeat(0.0, len(self.specifieds_symbols))))
        d.update(self.specifieds)
        return d

    @property
    def times(self):
        """An array-like object, containing time values over which the
        equations of motion are integrated, numerically.

        The object should be in a format which the integration module to be
        used can accept.
        """
        return self._times

    @times.setter
    def times(self, new_times):
        self._times = np.asarray(new_times)
        self._check_times(self._times)

    def _check_times(self, times):

        if len(times.shape) == 0:
            raise TypeError("Times supplied should be in an array_like format.")

        if not np.all(times >= 0):
            raise ValueError("Times supplied must have positive values.")

        if not np.all(np.diff(times) >= 0):
            raise ValueError("Times supplied should be in an ascending order.")

        return True

    @property
    def ode_solver(self):
        """A function that performs forward integration. It must have the
        same signature as scipy.integrate.odeint, which is::

            x_history = ode_solver(f, x0, t, args=(args,))

        where f is a function f(x, t, args), x0 are the initial conditions,
        x_history is the state time history, x is the state, t is the time,
        and args is a keyword argument takes arguments that are then passed
        to f. The default solver is odeint.

        """
        return self._ode_solver

    @ode_solver.setter
    def ode_solver(self, ode_solver):
        if not hasattr(ode_solver, '__call__'):
            msg = "``ode_solver`` ({}) is not a function."
            raise ValueError(msg.format(ode_solver))
        self._ode_solver = ode_solver

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
        d = dict(zip(self.states,
                     repeat(0.0, len(self.states))))
        d.update(self.initial_conditions)
        return d

    @property
    def _initial_conditions_array(self):
        x0_dict = self._initial_conditions_padded_with_defaults()
        return np.array([x0_dict[xi] for xi in self.states])

    @property
    def evaluate_ode_function(self):
        """A function generated by ``generate_ode_function`` that computes
        the state derivatives:

            x' = evaluate_ode_function(x, t, *args)

        This function is used by the ``ode_solver``.

        """
        return self._evaluate_ode_function

    def _args_for_gen_ode_func(self):
        """Returns a tuple of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generate_ode_function``.

        """

        args = (self.eom_method.forcing,
                self.coordinates,
                self.speeds,
                self.constants_symbols)

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

        kwargs = {
            'mass_matrix': self.eom_method.mass_matrix,
            'coordinate_derivatives': kin_diff_rhs,
            'specifieds': specifieds,
        }

        return kwargs

    def _generate_constraint_functions(self):
        # TODO : Support explicit instance of time in the equations.

        all_constraints = []

        if self.eom_method._f_h:
            self._eval_holonomic = sm.lambdify(
                (self.states, self.constants_symbols),
                self.eom_method._f_h, cse=True)
            all_constraints += self.eom_method._f_h[:]

        if self.eom_method._k_nh:
            # rebuild the nonholonomic constraints
            # TODO : KanesMethod and _Method should store the original
            # constraints passed by the user. Fix in sympy.physics.mechanics.
            nh = (self.eom_method._k_nh*self.eom_method.u +
                  self.eom_method._f_nh)
            self._eval_nonholonomic = sm.lambdify(
                (self.states, self.constants_symbols),
                nh, cse=True)
            all_constraints += nh[:]

        if all_constraints:
            all_constraints = sm.Matrix(all_constraints)
            self._eval_constraints = sm.lambdify(
                (self.states, self.constants_symbols),
                all_constraints, cse=True)
            self._all_constraints = all_constraints

    def evaluate_holonomic(self, x=None):
        """Returns the value of the holonomic constraints given the system
        state.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> system.generate_ode_function()
            >>> help(system.evaluate_ode_function)

        """
        if not self.eom_method._f_h:
            raise ValueError('This system has no holonomic constraints.')
        if not hasattr(self, '_eval_holonomic'):
            self._generate_constraint_functions()
        if x is None:
            x = self._initial_conditions_array
        else:
            x = np.asarray(x)
        p = self._constants_array
        if len(x.shape) == 2:
            x = x.T
            return self._eval_holonomic(x, p).squeeze().T
        else:
            return self._eval_holonomic(x, p).squeeze()

    def evaluate_nonholonomic(self, x=None):
        """Returns the value of the nonholonomic constraints given the system
        state.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> system.generate_ode_function()
            >>> help(system.evaluate_ode_function)

        """
        if not self.eom_method._k_nh:
            msg = 'This system has no nonholonomic constraints.'
            raise ValueError(msg)
        if not hasattr(self, '_eval_nonholonomic'):
            self._generate_constraint_functions()
        if x is None:
            x = self._initial_conditions_array
        else:
            x = np.asarray(x)
        p = self._constants_array
        if len(x.shape) == 2:
            x = x.T
            return self._eval_nonholonomic(x, p).squeeze().T
        else:
            return self._eval_nonholonomic(x, p).squeeze()

    def evaluate_constraints(self, x=None):
        """Returns the values of the holonomic and nonholonomic constraints at
        the initial condition or, alternatively, for the provided state vector.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State vector of n states or a series of m state vectors.

        Returns
        =======
        ndarray, shape(o,) or shape(m, o)
            Constraint vector of o constraints or a series of m constraint
            vectors.

        Notes
        =====

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> system.generate_ode_function()
            >>> help(system.evaluate_ode_function)

        """
        if not self.eom_method._f_h and not self.eom_method._k_nh:
            msg = 'This system does not have constraints.'
            raise ValueError(msg)
        if not hasattr(self, '_eval_constraints'):
            self._generate_constraint_functions()
        if x is None:
            x = self._initial_conditions_array
        else:
            x = np.asarray(x)
        p = self._constants_array
        if len(x.shape) == 2:
            x = x.T
            return self._eval_constraints(x, p).squeeze().T
        else:
            return self._eval_constraints(x, p).squeeze()

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
            Extra keyword arguments that are passed to ``scipy.optimize.root``.

        """
        # TODO : The nonholonomic constraints can be solved analytically, root
        # is only required for the coordinates.

        if not self.eom_method._f_h and not self.eom_method._k_nh:
            msg = ('This system does not have constraints, set all initial '
                   'conditions yourself.')
            raise ValueError(msg)
        else:
            num_holo = len(self.eom_method._f_h)
            num_nonh = self.eom_method._k_nh.shape[0]

        # TODO : These should be publicly accessible on KanesMethod.
        if dep_vars is None:
            dep_vars = self.eom_method._qdep[:] + self.eom_method._udep[:]

        # TODO : Would be nice to check if the dependent variables are present
        # in the constraints and that the right number of coordinates and
        # speeds are each supplied.
        if len(dep_vars) != num_holo + num_nonh:
            msg = (f'You must supply {num_holo} dependent coordinates and '
                   f'{num_nonh} dependent speeds.')
            raise ValueError(msg)

        x = self._initial_conditions_array
        p = self._constants_array
        dep_guess = [self.initial_conditions[xi] for xi in dep_vars]
        dep_idxs = [self.states.index(xi) for xi in dep_vars]

        if use_jac:
            jac = self._all_constraints.jacobian(dep_vars)
            eval_jac = sm.lambdify((self.states, self.constants_symbols), jac,
                                   cse=True)

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                return self.evaluate_constraints(x=x), eval_jac(x, p)

            fprime = True
        else:

            def eval_f(x_dep, p):
                x[dep_idxs] = x_dep
                return self.evaluate_constraints(x=x)

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
        """Calls ``pydy.codegen.ode_function_generators.generate_ode_function``
        with the appropriate arguments, and sets the ``evaluate_ode_function``
        attribute to the resulting function.

        Parameters
        ----------
        kwargs
            All other kwargs are passed onto
            ``pydy.codegen.ode_function_generators.generate_ode_function()``.
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
        self._check_times(self.times)

        if len(self.times) < 2:
            raise ValueError('The times vector must be at least length 2.')

        if self.evaluate_ode_function is None:
            self.generate_ode_function()

        init_conds_dict = self._initial_conditions_padded_with_defaults()
        initial_conditions_in_proper_order = [
            init_conds_dict[k] for k in self.states]

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

        return initial_conditions_in_proper_order, args

    def evaluate_ode(self, x=None, t=None):
        """Returns the right hand side of the differential equations. The
        default is to evaluate at the set initial_conditions at the first time
        value. Pass in optional arguments to override using the initial state
        and time.

        Parameters
        ==========
        x : array_like, shape(n,) or shape(m, n), optional
            State values at time t.
        t : float or array_like, shape(m,), optional
            Time or m time values.

        Returns
        =======
        x_dot : ndarray, shape(n,) or shape(m, n)
           Time derivative of the states at time t.

        Notes
        =====

        This method is present for convenience, it is not designed to be used
        where performance matters, use :py:meth:`System.evaluate_ode_function`
        directly when performance is needed.

        To see the order of the state values use::

            >>> system = System(...)
            >>> system.states

        or::

            >>> system.generate_ode_function()
            >>> help(system.evaluate_ode_function)

        """
        x_default, args = self._prep_for_evaluate()

        if x is None:
            x = x_default
        x = np.asarray(x)

        if t is None:
            if len(x.shape) == 1:
                t = self.times[0]
            else:
                t = self.times

        if len(x.shape) == 1 and not isinstance(t, float):
            raise ValueError('Time must be a float.')
        elif len(x.shape) == 1 and isinstance(t, float):
            return self.evaluate_ode_function(x, t, *args)

        # NOTE : I tried to make use of numpy.vectorize but it is not possible
        # due to args not being necessarily being comprised of arrays.
        if len(x.shape) == 2:
            if x.shape[0] != len(t):
                raise ValueError('x trajectory must have same length as t.')
            xdot = np.zeros_like(x)
            for i, (ti, xi) in enumerate(zip(t, x)):
                xdot[i, :] = self.evaluate_ode_function(xi, ti, *args)
            return xdot


    def integrate(self, **solver_kwargs):
        """Integrates the equations ``evaluate_ode_function()`` using
        ``ode_solver``.

        It is necessary to have first generated an ode function. If you have
        not done so, we do so automatically by invoking
        ``generate_ode_function()``. However, if you want to customize how this
        function is generated (e.g., change the generator to cython), you can
        call ``generate_ode_function()`` on your own (before calling
        ``integrate()``).

        Parameters
        ----------
        **solver_kwargs
            Optional arguments that are passed on to the ``ode_solver``.

        Returns
        -------
        x_history : np.array, shape(num_integrator_time_steps, 2)
            The trajectory of states (coordinates and speeds) through the
            requested time interval. num_integrator_time_steps is either
            len(times) if len(times) > 2, or is determined by the
            ``ode_solver``.

        """
        x0, args = self._prep_for_evaluate()

        # NOTE : User cannot pass in args, System handles that.
        solver_kwargs.pop('args', None)

        x_history = self.ode_solver(
            self.evaluate_ode_function,
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

    def _Kane_undefined_dynamicsymbols(self):
        """Similar to ``_find_dynamicsymbols()``, except that it checks all
        syms used in the system. Code is copied from ``linearize()``.

        TODO temporarily here until KanesMethod and Lagranges method have an
        interface for obtaining these quantities.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        functions_of_time = set()
        for expr in from_eoms:
            functions_of_time = functions_of_time.union(
                find_dynamicsymbols(expr))
        return functions_of_time.difference(from_sym_lists)

    def _Kane_constant_symbols(self):
        """Similar to ``_find_othersymbols()``, except it checks all syms used in
        the system.

        Remove the time symbol.

        TODO temporary.

        """
        from_eoms, from_sym_lists = self._Kane_inlist_insyms()
        unique_symbols = set()
        for expr in from_eoms:
            unique_symbols = unique_symbols.union(expr.free_symbols)
        constants = unique_symbols
        constants.remove(dynamicsymbols._t)
        return constants
