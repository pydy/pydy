
from .codegen.code import generate_ode_function
from scipy.integrate import odeint

class System(object):
    """Manages the simulation (integration) of a system whose equations are
    given by KanesMethod or LagrangesMethod.

    All attributes can be set directly. With the exception of `method`, the
    attributes can also be set via keyword arguments to the constructor.

    Attributes
    ----------
    method : sympy.physics.mechanics.KanesMethod or
                sympy.physics.mechanics.LagrangesMethod
        The method used to generate the equations of motion.
    constants : dict, optional (default: all 1.0)
        Numerical values for the constants in the problem. Keys are the symbols
        for the constants, and values are floats. Constants that are not
        specified in this dict are given a default value of 1.0.
    specified_symbols : iterable of symbols, optional
        The symbols for the quantities whose numerical values you want to
        specify. The order of these must match the order of the
        `specified_values`. If not provided, we get all the specified symbols
        from the `method`, and set their numerical value to 0. You can also
        provide a subset of symbols for this attribute, in which case the
        remaining symbols still have their numerical value set to 0.
    specified_values : iterable of floats or functions, optional
        The numerical values of the `specified_symbols`, or functions (which
        take the same arguments as f, above) that can generate the numerical
        values. If not provided, we get all the specified symbols
        from the `method`, and set their numerical value to 0. You can also
        provide a subset of symbols for this attribute, in which case the
        remaining symbols still have their numerical value set to 0.
    code_gen_backend : str, optional (default: 'lambdify')
        The backend used to generate the ode function.
        See the documentation of pydy.codegen.code.generate_ode_function`.
    ode_solver : function, optional (default: scipy.integrate.odeint)
        A function that performs forward integration. It must have the same
        signature as odeint, which is::
       
            x_history = ode_solver(f, x0, t, args=(args,))

        where f is a function f(x, t), x0 are the initial conditions, x_history
        is the history, and args is a keyword argument takes arguments that are
        then passed to f. TODO
    initial_conditions : dict, optional (default: all zero)
        Initial conditions for all coordinates and speeds. Keys are the symbols
        for the coordinates and speeds, and values are floats. Coordinates or
        speeds that are not specified in this dict are given a default value of
        zero.

    """
    def __init__(self, method, **kwargs):
        self.method = method
        for k, v in kwargs:
            setattr(self, k, v)
        if self.constants != None:
            self.constants = self._find_constants()
        if self.specified_symbols != None and self.specified_values != None:
            self._check_specified()
            # TODO
            pass
        if self.code_gen_backend != None:
            self.code_gen_backend = 'lambdify'
        if ode_solver != None:
            self.ode_solver = odeint

    def _find_constants(self):
        othersymbols = method._find_othersymbols()
        return zip(othersymbols, len(othersymbols) * [1.0])

    def _find_specified_symbols(self):
        specified_symbols = method._find_dynamicsymbols()
        return zip

    def generate_ode_function(self):
        """Calls `pydy. TODO
        
        """
        self.rhs = generate_ode_function(
                self.method.mass_matrix_full, self.method.forcing_full,
                self.constants.keys(),
                self.method._q, self.method._u)



#for k, v in kwargs:
#    setattr(self, k, v)
#sys.generate_ode_function(backend=...)
#sys.code_gen_backend = 
#sys.ode_solver = odeint
#sys.specified_symbols =
#sys.specified_values =
#sys.initial_conditions =
#sys.constants = 
#
#class Sys:
#    self.initial_conditions = zero(
#    def integrate(self, times):
#        # make rhs
#        sys.ode_solver(self.rhs, self.initial_conditions, times,
#                args=(self.con)
#    

