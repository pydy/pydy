
from .codegen.code import generate_ode_function
from scipy.integrate import odeint

class System(object):
    """Manages the simulation (integration) of a system whose equations are
    given by KanesMethod or LagrangesMethod.

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
        times = [0, 5.0]
        sys.integrate(times)

    In this case, we use defaults for the numerical values of the constants,
    specified quantities, initial conditions, etc. You probably won't like
    these defaults. You can also specify such values via constructor keyword
    arguments or via the attributes::

        sys = System(km, initial_conditions={symbol('q1'): 0.5})
        sys.constants = {symbol('m'): 5.0}
        sys.integrate(times)

    In this case, we generate the numerical ode function for you behind the
    scenes. If you want to customize how this function is generated, you must
    call `generate_ode_function` on your own::

        sys = System(KM)
        sys.generate_ode_function(generator='cython')
        sys.integrate(times)

    """

    default_ode_solver = odeint

    def __init__(self, eom_method, **kwargs):
        """See the class's attributes for a description of the arguments to
        this constructor.

        The parameters to this constructor are all attributes of the System.
        Actually, they are properties. With the exception of `eom_method`,
        these attributes can be modified directly at any future point.

        Parameters
        ----------
        eom_method : sympy.physics.mechanics.KanesMethod or
                sympy.physics.mechanics.LagrangesMethod
            If using `KanesMethod`, you must have called
            `KanesMethod.kanes_equations()` *before* constructing this
            `System`. Similarly, if using `LagrangesMethod`, you must have
            called `LagrangesMethod.form_lagranges_equations()`.
        constants : dict, optional (default: all 1.0)
        specifieds : dict, optional (default: all 0)
        ode_solver : function, optional (default: scipy.integrate.odeint)
        initial_conditions : dict, optional (default: all zero)

        """
        self._eom_method = eom_method

        # This will invoke the setters and do the necessary error-checking, and
        # fill in defaults.
        for k, v in kwargs:
            setattr(self, k, v)

        # For all attributes not provided in kwargs, set attribute to default.
        if self.constants == None: 
            # This will invoke the setting of defaults.
            self.constants = dict()

        if self.specifieds == None:
            # This will invoke the setting of defaults.
            self.specifieds = dict()

        if self.ode_solver == None:
            self.ode_solver = self.default_ode_solver

        if self.initial_conditions == None:
            # This will invoke the setting of defaults.
            self.initial_conditions = dict()

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
        """ This is either a sympy.physics.mechanics.KanesMethod or
        sympy.physics.mechanics.LagrangesMethod. The method used to generate
        the equations of motion. Read-only.

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

    @specifieds.setter
    def constants(self, constants):
        self._check_constants(constants)
        self._constants = self._default_constants()
        self._constants.update(constants)

    def _system_constants_symbols(self):
        """Wrapper."""
        return method._find_othersymbols()

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
            sys.specifieds = {(a, b, c): lambda x, t: k * x}

        """
        return self._specifieds

    @specifieds.setter
    def specifieds(self, specifieds):
        self._check_specifieds(specifieds)
        self._specifieds = self._default_specifieds()
        self._specifieds.update(specifieds)

    def _system_specifieds_symbols(self):
        """Wrapper."""
        return method._find_dynamicsymbols()

    def _default_specifieds(self):
        symbols = self._system_specifieds()
        return dict(zip(symbols, len(symbols) * [0.0]))
        
    def _assert_is_specified_symbol(symbol, all_symbols):
        if symbol not in all_symbols:
            raise Exception("Symbol {} is not a 'specified' symbol.".format(k))

    def _check_specifieds(self, specifieds):
        symbols = self._system_specifieds_symbols()

        for k, v in specifieds.items():

            # The symbols must be specifieds.
            if isinstance(k, tuple):
                for ki in k:
                    self._assert_is_specified_symbol(ki, symbols)
            else:
                self._assert_is_specified_symbol(k, symbols)

            # Length of the numerical values is same as number of symbols.
            if hasattr(v, '__call__'):
                # TODO how do we reliably check the output size of the 
                # function? Do we call the function with test inputs? That is
                # not reliable... My vote is that we just don't do the check.
                pass
            else:
                if len(k) != len(v):
                    raise ValueError("{} is not the same length as {}".format(
                        k, v))

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

    @initial_conditions.setter(self, initial_conditions):
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
                self.constants.keys(),
                self.coordinates, self.speeds,
                # kwargs:
                specified=self.specifieds.keys(),
                generator=generator,
                **kwargs,
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
        if self.evaluate_ode_function == None:
            self.generate_ode_function()
        # TODO double-check this.
        initial_conditions_in_proper_order = \
                [self.initial_conditions[k] for k in self.states]
        return self.ode_solver(
                self.evaluate_ode_function,
                initial_conditions_in_proper_order,
                times,
                args={
                    'constants': self.constants.values(),
                    'specified': self.specifieds.values()
                    }
                )