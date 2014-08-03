
from .codegen.code import generate_ode_function
from scipy.integrate import odeint
from sympy import symbols

class System(object):
    """Manages the simulation (integration) of a system whose equations are
    given by KanesMethod.

    Many of the attributes are also properties, and can be directly modified.

    The usage of this class is to:

        1. specify your options either via the constructor or via the
           attributes.
        2. optionally, call `generate_ode_function` if you want to customize
           how this function is generated.
        3. call `integrate` to simulate your system.

    Examples
    --------
    The simplest usage of this class is as follows::

        km = KanesMethod(...)
        km.kanes_equations(force_list, body_list)
        sys = System(km)
        times = np.linspace(0, 5, 100)
        sys.integrate(times)

    In this case, we use defaults for the numerical values of the constants,
    specified quantities, initial conditions, etc. You probably won't like
    these defaults. You can also specify such values via constructor keyword
    arguments or via the attributes::

        sys = System(km, initial_conditions={symbol('q1'): 0.5})
        sys.constants = {symbol('m'): 5.0}
        sys.integrate(times)

    In this case, the System generates the numerical ode function for you
    behind the scenes. If you want to customize how this function is generated,
    you must call `generate_ode_function` on your own::

        sys = System(KM)
        sys.generate_ode_function(generator='cython')
        sys.integrate(times)

    """
    def __init__(self, eom_method, constants=dict(), specifieds=dict(),
            ode_solver=odeint, initial_conditions=dict()):
        """See the class's attributes for a description of the arguments to
        this constructor.

        The parameters to this constructor are all attributes of the System.
        Actually, they are properties. With the exception of `eom_method`,
        these attributes can be modified directly at any future point.

        Parameters
        ----------
        eom_method : sympy.physics.mechanics.KanesMethod
            You must have called `KanesMethod.kanes_equations()` *before*
            constructing this `System`.
        constants : dict, optional (default: all 1.0)
        specifieds : dict, optional (default: all 0)
        ode_solver : function, optional (default: scipy.integrate.odeint)
        initial_conditions : dict, optional (default: all zero)

        """
        self._eom_method = eom_method
        self.constants = constants
        self.specifieds = specifieds
        self.ode_solver = ode_solver
        self.initial_conditions = initial_conditions
        self._evaluate_ode_function = None

    @property
    def coordinates(self):
        return self.method._q

    @property
    def sppeds(self):
        return self.method._u

    @property
    def states(self):
        """These are in the same order as used in integration (as passed into
        evaluate_ode_function).

        """
        return self.method._q + self.method._u

    @property
    def eom_method(self):
        """ This is a sympy.physics.mechanics.KanesMethod. The method used to
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
        self._constants = self._default_constants()
        self._constants.update(constants)

    def _system_constants_symbols(self):
        """Wrapper."""
        # TODO is it expensive to make repeated calls to this?
        return self._Kane_constant_symbols()

    def _default_constants(self):
        symbols = self._system_constants_symbols()
        return dict(zip(symbols, len(symbols) * [1.0]))

    def _check_constants(self, constants):
        symbols = self._system_constants_symbols()
        for k in constants.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a constant.".format(k))

    @property
    def specifieds(self):
        """A dict that provides numerical values for the specified quantities
        in the problem (all dynamicsymbols that are not given by the `method`).
        Keys are the symbols for the specified quantities, or a tuple of
        symbols, and values are the floats, arrays of floats, or functions that
        generate the values. If a dictionary value is a function, it must have
        the same signature as `f(x, t)`, the ode right-hand-side function (see
        the documentation for the `ode_solver` attribute). You needn't provide
        values for all specified symbols. Those for which you do not give a
        value will default to 0.0.

        Examples
        --------
        Keys can be individual symbols, or a tuple of symbols. Length of a
        value must match the length of the corresponding key. Values can be
        functions that return iterables::

            sys = System(km)
            sys.specifieds = {(a, b, c): np.ones(3), d: lambda x, t: -3 * x[0]}
            sys.specifieds = {(a, b, c): lambda x, t: np.ones(3)}

        """
        return self._specifieds

    @specifieds.setter
    def specifieds(self, specifieds):
        self._check_specifieds(specifieds)
        self._specifieds = specifieds
        self._fill_in_default_specifieds()

    def _system_specifieds_symbols(self):
        """Wrapper."""
        # TODO is it expensive to make repeated calls to this?
        # TODO eventually use a method in the KanesMethod class.
        return self._Kane_undefined_dynamicsymbols()

    def _default_specifieds(self):
        symbols = self._system_specifieds_symbols()
        return dict(zip(symbols, len(symbols) * [0.0]))
        
    def _assert_is_specified_symbol(self, symbol, all_symbols):
        if symbol not in all_symbols:
            raise ValueError("Symbol {} is not a 'specified' symbol.".format(k))

    def _assert_symbol_appears_multiple_times(self, symbol, specifieds):
        if symbol in specifieds:
            raise ValueError("Symbol {} appears in {} more than once.".format(
                symbol, specifieds))

    def _check_specifieds(self, specifieds):
        symbols = self._system_specifieds_symbols()

        symbols_so_far = list()

        for k, v in specifieds.items():

            # The symbols must be specifieds.
            if isinstance(k, tuple):
                for ki in k:
                    self._assert_is_specified_symbol(ki, symbols)
            else:
                self._assert_is_specified_symbol(k, symbols)

            # Length of the numerical values is same as number of symbols.
            if hasattr(v, '__call__'):
                # How do we reliably check the output size of the 
                # function? Do we call the function with test inputs? That is
                # not reliable... My vote is that we just don't do the check.
                pass
            else:
                if len(k) != len(v):
                    raise ValueError("{} is not the same length as {}".format(
                        k, v))

            # Each specified symbol can appear only once.
            if isinstance(k, tuple):
                for ki in k:
                    self._assert_symbol_appears_multiple_times(ki, specifieds)
                    symbols_so_far.append(ki)
            else:
                self._assert_symbol_appears_multiple_times(ki, specifieds)
                symbols_so_far.append(k)

    def _symbol_is_in_specifieds_dict(self, symbol, specifieds_dict):
        for k in specifieds_dict.keys():
            if symbol == k:
                return True
            elif isinstance(k, tuple):
                return symbol in k
            else:
                raise ValueError("Unexpected key {}.".format(k))
        return False

    def _fill_in_default_specifieds(self):
        if self.specifieds == None:
            self.specifieds = dict()
        for symbol in self._system_specifieds_symbols():
            if not self._symbol_is_in_specifieds_dict(symbol, self.specifieds):
                self._specifieds[symbol] = 0.0

    @property
    def ode_solver(self):
        """A function that performs forward integration. It must have the same
        signature as odeint, which is::
       
            x_history = ode_solver(f, x0, t, args=(args,))

        where f is a function f(x, t), x0 are the initial conditions, x_history
        is the history, x is the state, t is the time, and args is a keyword
        argument takes arguments that are then passed to f. The default is
        odeint.

        """
        return self._ode_solver

    @ode_solver.setter
    def ode_solver(self, ode_solver):
        if not hasattr(ode_solver, '__call__'):
            raise ValueError(
                    "`ode_solver` ({}) is not a function.".format(ode_solver))
        self._ode_solver = ode_solver

    @property
    def initial_conditions(self):
        """ Initial conditions for all coordinates and speeds. Keys are the
        symbols for the coordinates and speeds, and values are floats.
        Coordinates or speeds that are not specified in this dict are given a
        default value of 0.0.

        """
        return self._initial_conditions

    @initial_conditions.setter
    def initial_conditions(self, initial_conditions):
        self._check_initial_conditions(initial_conditions)
        self._initial_conditions = self._default_initial_conditions()
        self._initial_conditions.update(initial_conditions)

    def _default_initial_conditions(self):
        symbols = self.states
        return dict(zip(symbols, len(symbols) * [0.0]))

    def _check_initial_conditions(self, constants):
        symbols = self.states
        for k in constants.keys():
            if k not in symbols:
                raise ValueError("Symbol {} is not a state.".format(k))

    @property
    def evaluate_ode_function(self):
        """A function generated by `generate_ode_function` that computes the
        state derivatives:
        
            x' = evaluate_ode_function(x, t, args=(...))

        This function is used by the `ode_solver`.

        """
        return self._evaluate_ode_function

    def generate_ode_function(self, generator='lambdify', **kwargs):
        """Calls `pydy.codegen.code.generate_ode_function` with the appropriate
        arguments, and sets the `evaluate_ode_function` attribute to the
        resulting function.

        Parameters
        ----------
        generator : str, optional (default: 'lambdify')
            See documentation for `pydy.codegen.code.generate_ode_function`
        kwargs
            All other kwargs are passed onto
            `pydy.codegen.code.generate_ode_function`. Don't specify the
            `specified` kwarg though; this class takes care of those.

        Returns
        -------
        evaluate_ode_function : function
            A function which evaluates the derivaties of the states.
        
        """
        self._evaluate_ode_function = generate_ode_function(
                # args:
                self.method.mass_matrix_full, self.method.forcing_full,
                self._system_constants_symbols(),
                self.coordinates, self.speeds,
                # kwargs:
                specified=self._system_specifieds_symbols(),
                generator=generator,
                **kwargs
                )
        return self.evaluate_ode_function

    def integrate(self, times):
        """Integrates the equations `evaluate_ode_function` using `ode_solver`.

        It is necessary to have first generated an ode function. If you have
        not done so, we do so automatically by invoking
        `generate_ode_function`. However, if you want to customize how this
        function is generated (e.g., change the generator to cython), you can
        call `generate_ode_function` on your own (before calling `integrate`).

        Returns
        -------
        x_history : np.array, shape(num_integrator_time_steps, 2)
            The trajectory of states (coordinates and speeds) through the
            requested time interval. num_integrator_time_steps is either
            len(times) if len(times) > 2, or is determined by the `ode_solver`.

        """
        # Users might have changed these properties by directly accessing the
        # dict, without using the setter. Before we integrate, make sure they
        # did not muck up these dicts.
        _check_constants(self.constants)
        _check_specifieds(self.specifieds)
        _check_initial_conditions(self.initial_conditions)

        if self.evaluate_ode_function == None:
            self.generate_ode_function()
        initial_conditions_in_proper_order = \
                [self.initial_conditions[k] for k in self.states]
        return self.ode_solver(
                self.evaluate_ode_function,
                initial_conditions_in_proper_order,
                times,
                args={
                    'constants': self.constants,
                    'specified': self.specifieds,
                    }
                )

    def _Kane_inlist_insyms(self):
        """TODO temporary."""
        uaux = self._eom_method._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds & derivatives which are equal to zero
        subdict = dict(
                list(zip(uaux + uauxdot, [0] * (len(uaux) + len(uauxdot)))))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        insyms = set(self._eom_method._q + self._eom_method._qdot +
                self._eom_method._u + self._eom_method._udot + uaux + uauxdot)
        inlist = self._eom_method._f_d.subs(subdict)
        return inlist, insyms

    def _Kane_undefined_dynamicsymbols(self):
        """Similar to `_find_dynamicsymbols()`, except that it checks all syms
        used in the system. Code is copied from `linearize()`.

        TODO temporarily here until KanesMethod and Lagranges method have an
        interface for obtaining these quantities.

        """
        return list(self._eom_method._find_dynamicsymbols(
            *self._Kane_inlist_insyms()))

    def _Kane_constant_symbols(self):
        """Similar to `_find_othersymbols()`, except it checks all syms used in
        the system.

        Remove the time symbol.

        TODO temporary.

        """
        constants = list(self._eom_method._find_othersymbols(
            *self._Kane_inlist_insyms()))
        constants.remove(symbols('t'))
        return constants

   
