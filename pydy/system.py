"""The System class manages the simulation (integration) of a system whose
equations are given by KanesMethod.

Many of the attributes are also properties, and can be directly modified.

Here is the procedure for using this class.

    1. specify your options either via the constructor or via the
       attributes.
    2. optionally, call ``generate_ode_function()`` if you want to customize
       how the ODE function is generated.
    3. call ``integrate()`` to simulate your system.

The simplest usage of this class is as follows. First, we need a
KanesMethod object on which we have already invoked ``kanes_equations()``::

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
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics.functions import find_dynamicsymbols
from scipy.integrate import odeint

from .codegen.ode_function_generators import generate_ode_function
from .utils import PyDyFutureWarning

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
    specifieds : dict, optional (default: all 0)
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
            self._times = []
        else:
            self._times = times

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

         Here are examples for (2):

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
    def evaluate_ode_function(self):
        """A function generated by ``generate_ode_function`` that computes
        the state derivatives:

            x' = evaluate_ode_function(x, t, args=(...))

        This function is used by the ``ode_solver``.

        """
        return self._evaluate_ode_function

    def _args_for_gen_ode_func(self):
        """Returns a tuple of arguments in the form required by
        ``pydy.codegen.ode_function_generators.generate_ode_function``.

        """

        args = (self.eom_method.forcing_full,
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

        kwargs = {'mass_matrix': self.eom_method.mass_matrix_full,
                  'specifieds': specifieds}

        return kwargs

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

        """

        if 'specified' in kwargs:
            kwargs.pop('specified')
            print("User supplied 'specified' kwarg was disregarded.")

        if 'specifieds' in kwargs:
            kwargs.pop('specifieds')
            print("User supplied 'specifieds' kwarg was disregarded.")

        kwargs.update(self._kwargs_for_gen_ode_func())

        self._evaluate_ode_function = generate_ode_function(
            *self._args_for_gen_ode_func(),
            **kwargs)

        return self.evaluate_ode_function

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
        # Users might have changed these properties by directly accessing the
        # dict, without using the setter. Before we integrate, make sure they
        # did not muck up these dicts.
        self._check_constants(self.constants)
        self._check_specifieds(self.specifieds)
        self._check_initial_conditions(self.initial_conditions)
        self._check_times(self.times)

        if self.evaluate_ode_function is None:
            self.generate_ode_function()

        init_conds_dict = self._initial_conditions_padded_with_defaults()
        initial_conditions_in_proper_order = \
            [init_conds_dict[k] for k in self.states]

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

        x_history = self.ode_solver(
            self.evaluate_ode_function,
            initial_conditions_in_proper_order,
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
