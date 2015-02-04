#!/usr/bin/env python

import numpy as np
import numpy.linalg
import scipy.linalg
import sympy as sm
import sympy.physics.mechanics as me
from sympy.core.function import UndefinedFunction
Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')
if theano:
    from sympy.printing.theanocode import theano_function

from .cython_code import CythonMatrixGenerator


class ODEFunctionGenerator(object):
    """This is an abstract base class for all of the generators. A subclass
    is expected to implement the methods necessary to evaluate the arrays
    needed to compute xdot for the three different system specification
    types."""

    _rhs_doc_template = \
"""\
Returns the derivatives of the states, i.e. numerically evaluates the right
hand side of the first order differential equation.

x' = f(x, t,{specified_call_sig} p)

Parameters
==========
x : ndarray, shape({num_states},)
    The state vector is ordered as such:
{state_list}
t : float
    The current time.{specified_explanation}
p : dictionary len({num_constants}) or ndarray shape({num_constants},)
    The dictionary maps the constant symbols to floats. The constants are,
    in order:
{constant_list}

Returns
=======
dx : ndarray, shape({num_states},)
    The derivative of the state vector.

"""

    _specifieds_doc_template = \
"""
r : dictionary; ndarray, shape({num_specified},); function

    There are three options for this dictionary. (1) is more flexible but
    (2) and (3) are much more efficient.

    (1) A dictionary that maps the specified functions of time to floats,
    ndarrays, or functions that produce ndarrays. The keys can be a single
    specified symbolic function of time or a tuple of symbols. The total
    number of symbols must be equal to {num_specified}. If the value is a
    function it must be of the form g(x, t), where x is the current state
    vector ndarray and t is the current time float and it must return an
    ndarray of the correct shape. For example::

      r = {{a: 1.0,
           (d, b) : np.array([1.0, 2.0]),
           (e, f) : lambda x, t: np.array(x[0], x[1]),
           c: lambda x, t: np.array(x[2])}}

    (2) A ndarray with the specified values in the correct order and of the
    correct shape.

    (3) A function that must be of the form g(x, t), where x is the current
    state vector and t is the current time and it must return an ndarray of
    the correct shape.

    The specified inputs are, in order:
{specified_list}\
"""

    @staticmethod
    def _deduce_system_type(**kwargs):
        """Based on the combination of arguments this returns which ODE
        description has been provided.

        full rhs
            x' = f(x, t, r, p)
        full mass matrix
            M(x, p) * x' = f(x, t, r, p)
        min mass matrix
            M(q, p) * u' = f(q, u, t, r, p)
            q' = g(q, u, t)

        """

        if kwargs.pop('coordinate_derivatives') is not None:
            system_type = 'min mass matrix'
        elif kwargs.pop('mass_matrix') is not None:
            system_type = 'full mass matrix'
        else:
            system_type = 'full rhs'

        return system_type

    def __init__(self, right_hand_side, coordinates, speeds, constants,
                 mass_matrix=None, coordinate_derivatives=None,
                 specifieds=None, linear_sys_solver='numpy',
                 constants_arg_type=None, specifieds_arg_type=None):
        """

        Parameters
        ==========
        right_hand_side : SymPy Matrix, shape(n, 1)
            A column vector containing the symbolic expressions for the
            right hand side of the ordinary differential equations. If the
            right hand side has been solved for symbolically then x' = f(x,
            t, r, p), if it hasn't then you must supply the accompanying
            mass matrix M(x) * x' = f(x, t, r, p).
        coordinates : sequence of SymPy Functions
            The generalized coordinates.
        speeds : sequence of SymPy Functions
            The generalized speeds.
        constants : sequence of SymPy Symbols
            All of the constants present in the right hand side and mass
            matrix.
        mass_matrix : sympy.Matrix, shape(n, n), optional
            This can be either the "full" mass matrix in M(x) * x' = f(x, t,
            r, p) or the minimal mass matrix in M(q) * u' = f(u, q, t, r,
            p).
        coordinate_derivatives : sympy.Matrix, shape(m, 1), optional
            If the "minimal" mass matrix is supplied, then this matrix
            represents the right hand side of q' = g(q, t, p).
        specifieds : sequence of SymPy Functions
            The specified exogenous inputs to the ODEs. These should be
            functions of time.
        linear_sys_solver : string or function
            Specify either `numpy` or `scipy` to use the linear solvers
            provided in each package or supply a function that solves a
            linear system Ax=b with the call signature x = solve(A, b). If
            you need to use custom kwargs for the scipy solver, pass in a
            lambda function that wraps the solver and sets them.
        constants_arg_type : string
            The generated function accepts two different types of arguments
            for the numerical values of the constants: an ndarray of the
            constants values in the correct order or a dictionary mapping
            the constants sysmbols to the numerical values. If None, this is
            determined inside of the generated function and can cause a
            significant slow down for performance critical code. If you know
            apriori what arg types you need to support choose either
            ``array`` or ``dictionary``. Note that ``array`` is faster than
            ``dictionary``.
        specifieds_arg_type : string
            The generated function accepts three different types of arguments
            for the numerical values of the specifieds: an ndarray of the
            specifieds values in the correct order, a function that
            generates the correct ndarray, a dictionary mapping the
            specifieds symbols or tuples of thereof to floats, ndarrays, or
            functions. If None, this is determined inside of the generated
            function and can cause a significant slow down for performance
            critical code. If you know apriori what arg types you want to
            support choose either ``array``, ``function``, or
            ``dictionary``. The speed of each, from fast to slow, are
            ``array``, ``function``, ``dictionary``.

        Notes
        =====
        The generated function still supports the pre-0.3.0 extra argument
        style, i.e. args = {'constants': ..., 'specified': ...}, but only if
        ``constants_arg_type`` and ``specifieds_arg_type`` are both set to
        None. This functionality is deprecated and will be removed in 0.4.0,
        so it's best to adjust your code to support the new argument types.
        See the docstring for the generated function for more info.

        """

        self.right_hand_side = right_hand_side
        self.coordinates = coordinates
        self.speeds = speeds
        self.constants = constants
        self.mass_matrix = mass_matrix
        self.coordinate_derivatives = coordinate_derivatives
        self.specifieds = specifieds
        self.linear_sys_solver = linear_sys_solver
        self.constants_arg_type = constants_arg_type
        self.specifieds_arg_type = specifieds_arg_type

        self.system_type = self._deduce_system_type(
            mass_matrix=mass_matrix,
            coordinate_derivatives=coordinate_derivatives)

        self.num_coordinates = len(coordinates)
        self.num_speeds = len(speeds)
        self.num_states = self.num_coordinates + self.num_speeds
        self.num_constants = len(constants)

        if self.specifieds is None:
            self.num_specifieds = 0
            self.specifieds_arg_type = None
        else:
            self.num_specifieds = len(specifieds)

        # These are pre-allocated storage for the numerical values used in
        # some of the rhs() evaluations.
        self._constants_values = np.empty(self.num_constants)
        self._specifieds_values = np.empty(self.num_specifieds)

        self._check_system_consitency()
        self._set_linear_sys_solver()

    def _set_linear_sys_solver(self):

        if isinstance(self.linear_sys_solver, type(lambda x: x)):
            self.solve_linear_system = self.linear_sys_solver
        elif self.linear_sys_solver == 'numpy':
            self.solve_linear_system = numpy.linalg.solve
        elif self.linear_sys_solver == 'scipy':
            self.solve_linear_system = scipy.linalg.solve
        else:
            msg = '{} is not a valid solver'
            raise ValueError(msg.format(self.linear_sys_solver))

    def _check_system_consitency(self):

        if self.system_type == 'min mass matrix':

            nr, nc = self.mass_matrix.shape
            assert self.num_speeds == nr == nc
            assert self.num_speeds == self.right_hand_side.shape[0]
            assert self.num_coordinates == self.coordinate_derivatives.shape[0]

        elif self.system_type == 'full mass matrix':

            nr, nc = self.mass_matrix.shape
            assert self.num_states == nr == nc
            assert self.num_states == self.right_hand_side.shape[0]
            assert self.coordinate_derivatives is None

        elif self.system_type == 'full rhs':

            assert self.num_states == self.right_hand_side.shape[0]
            assert self.mass_matrix is None
            assert self.coordinate_derivatives is None

    @staticmethod
    def list_syms(indent, syms):
        indentation = ' ' * indent
        lst = '- ' + ('\n' + indentation + '- ').join([str(s) for s in syms])
        return indentation + lst

    def _parse_old_style_extra_args(self, *args):
        """Returns the post-0.3.0 style args if the pre-0.3.0 style args are
        passed in. The pre-0.3.0 style args always have three args: (x, t,
        d) where d is is a dictionary which should always at least contain
        the key 'constants'. It may also contain a key 'specified'."""

        # DEPRECATED : Remove before 0.4.0 release.

        last_arg = args[-1]
        try:
            constants = last_arg['constants']
        # ValueError is needed for older NumPy versions.
        except (KeyError, IndexError, ValueError):
            return args
        else:
            new_args = list(args[:-1])  # gets x and t

            if self.specifieds is not None:
                new_args.append(last_arg['specified'])

            new_args.append(constants)

            return tuple(new_args)

    def _convert_constants_dict_to_array(self, p):
        """Returns an array of numerical values from the constants
        dictionary in the correct order."""

        # NOTE : It's unfortunate that this has to be run at every rhs eval,
        # because subsequent calls to rhs() doesn't require different
        # constants. I suppose you can sub out all the constants in the EoMs
        # before passing them into the generator. That would beg for the
        # capability to support self.constants=None to skip all of this
        # stuff in the rhs eval.
        for i, c in enumerate(self.constants):
            self._constants_values[i] = p[c]

        return self._constants_values

    def _parse_constants(self, *args):
        """Returns an ndarray containing the numerical values of the
        constants in the correct order. If the constants are already an
        array, that array is returned."""

        p = args[-1]
        try:
            p = self._convert_constants_dict_to_array(p)
        except IndexError:
            # p is an array so just return the args
            return args
        else:
            return args[:-1] + (p,)

    def _convert_specifieds_dict_to_array(self, x, t, r):

        for k, v in r.items():
            # TODO : Not sure if this is the best check here.
            if isinstance(type(k), UndefinedFunction):
                k = (k,)
            idx = [self.specifieds.index(symmy) for symmy in k]
            try:
                self._specifieds_values[idx] = v(x, t)
            except TypeError:  # not callable
                # If not callable, then it should be a float, ndarray,
                # or indexable.
                self._specifieds_values[idx] = v

        return self._specifieds_values

    def _parse_specifieds(self, x, t, r, p):

        if isinstance(r, dict):
            # NOTE : This function sets self._specifieds_values, so here we
            # return nothing.
            self._convert_specifieds_dict_to_array(x, t, r)
        else:
            # More efficient.
            try:
                self._specifieds_values[:] = r(x, t)
            except TypeError:  # not callable.
                # If not callable, then it should be a float or ndarray.
                self._specifieds_values[:] = r

        return x, t, self._specifieds_values, p

    def _parse_all_args(self, *args):
        """Returns args formatted for the post 0.3.0 generators using all of
        the parsers. This is the slowest method and is used by default if no
        information is provided by the user on which type of args will be
        passed in."""

        args = self._parse_old_style_extra_args(*args)

        args = self._parse_constants(*args)

        if self.specifieds is not None:
            args = self._parse_specifieds(*args)

        return args

    def _generate_rhs_docstring(self):

        template_values = {'num_states': self.num_states,
                           'state_list': self.list_syms(8, self.coordinates
                                                        + self.speeds),
                           'num_constants': self.num_constants,
                           'constant_list': self.list_syms(8, self.constants),
                           'specified_call_sig': '',
                           'specified_explanation': ''}

        if self.specifieds is not None:
            template_values['specified_call_sig'] = ' r,'
            specified_template_values = {
                'num_specified': self.num_specifieds,
                'specified_list': self.list_syms(8, self.specifieds)}
            template_values['specified_explanation'] = \
                self._specifieds_doc_template.format(
                    **specified_template_values)

        return self._rhs_doc_template.format(**template_values)

    def _create_rhs_function(self):
        """Returns a function in the form expected by scipy.integrate.odeint
        that computes the derivatives of the states."""

        # This god awful mess below exists because of the need to optimize
        # the speed of the rhs evaluation. We unfortunately support way too
        # many ways to pass in extra arguments to the generated rhs
        # function. The default behavior is to parse the arguments passed
        # into the rhs function which can add a lot of computational
        # overhead. So we allow the user to specify what type the extra args
        # should be for both the constants and the specifieds. The constants
        # can be None, 'array', or 'dictionary'. The specifieds can be None,
        # 'array', 'function', or 'dictionary'. Thus we have 12 permutations
        # of this "switch".

        p_arg_type = self.constants_arg_type
        r_arg_type = self.specifieds_arg_type

        def slice_x(x):
            q = x[:self.num_coordinates]
            u = x[self.num_coordinates:]
            return q, u

        if p_arg_type is None and r_arg_type is None:

            # This is the only rhs that will properly check for the
            # pre-0.3.0 rhs args for backwards compatibility.

            def rhs(*args):
                # args: x, t, p
                # or
                # args: x, t, r, p

                args = self._parse_all_args(*args)

                q, u = slice_x(args[0])

                xdot = self._base_rhs(q, u, *args[2:])

                return xdot

        elif p_arg_type == 'array' and r_arg_type is None:

            # This could be combined with:
            # elif p_arg_type == 'array' and r_arg_type == 'array':

            def rhs(*args):
                # args: x, t, p
                # or
                # args: x, t, r, p

                if self.specifieds is not None:
                    args = self._parse_specifieds(*args)

                q, u = slice_x(args[0])

                return self._base_rhs(q, u, *args[2:])

        elif p_arg_type == 'dictionary' and r_arg_type is None:

            # This could be combined with:
            # elif p_arg_type == 'dictionary' and r_arg_type == 'array':

            def rhs(*args):
                # args: x, t, p
                # or
                # args: x, t, r, p

                if self.specifieds is not None:
                    args = self._parse_specifieds(*args)

                p = self._convert_constants_dict_to_array(args[-1])

                q, u = slice_x(args[0])

                xdot = self._base_rhs(q, u, *(args[2:-1] + (p,)))

                return xdot

        # All of the cases below must have specifieds, so the number of args
        # is known. r_arg_type is forces to be None if self.specifieds is
        # None.

        elif p_arg_type is None and r_arg_type == 'array':

            def rhs(*args):
                # args: x, t, r, p

                args = self._parse_constants(*args)

                q, u = slice_x(args[0])

                return self._base_rhs(q, u, *args[2:])

        elif p_arg_type == 'array' and r_arg_type == 'array':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                return self._base_rhs(q, u, *args[2:])

        elif p_arg_type == 'dictionary' and r_arg_type == 'array':

            def rhs(*args):
                # args: x, t, r, p

                p = self._convert_constants_dict_to_array(args[-1])

                q, u = slice_x(args[0])

                return self._base_rhs(q, u, *(args[2:-1] + (p,)))

        elif p_arg_type is None and r_arg_type == 'dictionary':

            def rhs(*args):
                # args: x, t, r, p

                args = self._parse_constants(*args)

                q, u = slice_x(args[0])

                r = self._convert_specifieds_dict_to_array(*args[:3])

                return self._base_rhs(q, u, r, args[-1])

        elif p_arg_type == 'array' and r_arg_type == 'dictionary':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                r = self._convert_specifieds_dict_to_array(*args[:3])

                return self._base_rhs(q, u, r, args[-1])

        elif p_arg_type == 'dictionary' and r_arg_type == 'dictionary':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                p = self._convert_constants_dict_to_array(args[-1])

                r = self._convert_specifieds_dict_to_array(*args[:3])

                return self._base_rhs(q, u, r, p)

        elif p_arg_type is None and r_arg_type == 'function':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                args = self._parse_constants(*args)

                r = args[2](*args[:2])

                return self._base_rhs(q, u, r, args[-1])

        elif p_arg_type == 'array' and r_arg_type == 'function':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                r = args[2](*args[:2])

                return self._base_rhs(q, u, r, args[-1])

        elif p_arg_type == 'dictionary' and r_arg_type == 'function':

            def rhs(*args):
                # args: x, t, r, p

                q, u = slice_x(args[0])

                p = self._convert_constants_dict_to_array(args[-1])

                r = args[2](*args[:2])

                return self._base_rhs(q, u, r, p)

        # TODO : This should generate the appropriate docstring based on the
        # above if statements.
        rhs.__doc__ = self._generate_rhs_docstring()

        return rhs

    def _define_inputs(self):

        self.inputs = [self.coordinates, self.speeds, self.constants]

        if self.specifieds is not None:
            self.inputs.insert(2, self.specifieds)

    def _create_base_rhs_function(self):

        if self.system_type == 'full rhs':

            self._base_rhs = self.eval_arrays

        elif self.system_type == 'full mass matrix':

            def base_rhs(*args):

                M, F = self.eval_arrays(*args)
                return self.solve_linear_system(M, F)

            self._base_rhs = base_rhs

        elif self.system_type == 'min mass matrix':

            xdot = np.empty(self.num_states, dtype=float)

            def base_rhs(*args):
                M, F, qdot = self.eval_arrays(*args)
                if self.num_speeds == 1:
                    udot = F / M
                else:
                    udot = self.solve_linear_system(M, F)
                xdot[:self.num_coordinates] = qdot
                xdot[self.num_coordinates:] = udot
                return xdot

            self._base_rhs = base_rhs

    def generate(self):
        """Returns a function that evaluates the right hand side of the
        first order ordinary differential equations in one of two forms:

        x' = f(x, t, p)

        or

        x' = f(x, t, r, p)

        """

        if self.system_type == 'full rhs':
            self.generate_full_rhs_function()
        elif self.system_type == 'full mass matrix':
            self.generate_full_mass_matrix_function()
        elif self.system_type == 'min mass matrix':
            self.generate_min_mass_matrix_function()

        self._create_base_rhs_function()

        return self._create_rhs_function()


class CythonODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        if Cython is None:
            raise Exception('Cython must be installed to use this class.')
        else:
            super(CythonODEFunctionGenerator, self).__init__(*args, **kwargs)

    @staticmethod
    def _cythonize(outputs, inputs):
        return CythonMatrixGenerator(inputs, outputs).compile()

    def _set_eval_array(self, f):

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: f(q, u, p, *self.empties)
        else:
            self.eval_arrays = lambda q, u, r, p: f(q, u, r, p, *self.empties)

    def generate_full_rhs_function(self):

        self._define_inputs()
        outputs = [self.right_hand_side]

        self.empties = (np.empty(self.num_states, dtype=float),)

        self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_full_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        mass_matrix_result = np.empty(self.num_states ** 2, dtype=float)
        rhs_result = np.empty(self.num_states, dtype=float)

        self.empties = (mass_matrix_result, rhs_result)

        self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_min_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        mass_matrix_result = np.empty(self.num_speeds ** 2, dtype=float)
        rhs_result = np.empty(self.num_speeds, dtype=float)
        kin_diffs_result = np.empty(self.num_coordinates, dtype=float)
        self.empties = (mass_matrix_result, rhs_result, kin_diffs_result)

        self._set_eval_array(self._cythonize(outputs, self.inputs))


class LambdifyODEFunctionGenerator(ODEFunctionGenerator):

    def _lambdify(self, outputs):
        # TODO : We could forgo this substitution for speed purposes and
        # have lots of args for lambdify (like it used to be done) but there
        # may be some limitations on number of args.
        subs = {}
        vec_inputs = []
        if self.specifieds is None:
            def_vecs = ['q', 'u', 'p']
        else:
            def_vecs = ['q', 'u', 'r', 'p']

        for syms, vec_name in zip(self.inputs, def_vecs):
            v = sm.DeferredVector(vec_name)
            for i, sym in enumerate(syms):
                subs[sym] = v[i]
            vec_inputs.append(v)

        try:
            outputs = [me.msubs(output, subs) for output in outputs]
        except AttributeError:
            # msubs doesn't exist in SymPy < 0.7.6.
            outputs = [output.subs(subs) for output in outputs]

        modules = [{'ImmutableMatrix': np.array}, 'numpy']

        return sm.lambdify(vec_inputs, outputs, modules=modules)

    def generate_full_rhs_function(self):

        self._define_inputs()
        outputs = [self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: np.squeeze(f(q, u, p))
        else:
            self.eval_arrays = lambda q, u, r, p: np.squeeze(f(q, u, r, p))

    def generate_full_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                      f(q, u, p)])
        else:
            self.eval_arrays = lambda q, u, r, p: tuple([np.squeeze(o) for o
                                                         in f(q, u, r, p)])

    def generate_min_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                      f(q, u, p)])
        else:
            self.eval_arrays = lambda q, u, r, p: tuple([np.squeeze(o) for o
                                                         in f(q, u, r, p)])


class TheanoODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        if theano is None:
            raise Exception('Theano must be installed to use this class.')
        else:
            super(TheanoODEFunctionGenerator, self).__init__(*args, **kwargs)

    def _define_inputs(self):

        if self.specifieds is None:
            self.inputs = self.coordinates + self.speeds + self.constants
        else:
            self.inputs = (self.coordinates + self.speeds + self.specifieds
                           + self.constants)

    def _theanoize(self, outputs):

        self._define_inputs()

        f = theano_function(self.inputs, outputs, on_unused_input='ignore')

        # Theano will run faster if you trust the input. I'm not sure
        # what the implications of this are. See:
        # http://deeplearning.net/software/theano/tutorial/faq.html#faster-small-theano-function
        # Note that map(np.asarray, np.hstack(args)) is required if
        # trust_input is True. If it is False, then it will sanitize the
        # inputs. I'm not sure which one is faster.
        f.trust_input = True

        return f

    def generate_full_rhs_function(self):

        outputs = [self.right_hand_side]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return np.squeeze(f(*vals))

        self.eval_arrays = eval_arrays

    def generate_full_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return tuple([np.squeeze(o) for o in f(*vals)])

        self.eval_arrays = eval_arrays

    def generate_min_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return tuple([np.squeeze(o) for o in f(*vals)])

        self.eval_arrays = eval_arrays


def generate_ode_function(right_hand_side, coordinates, speeds, constants,
                          mass_matrix=None, coordinate_derivatives=None,
                          specifieds=None, linear_sys_solver='numpy',
                          constants_arg_type=None,
                          specifieds_arg_type=None,generator='lambdify'):
    """Returns a numerical function which can evaluate the right hand side
    of the first order ordinary differential equations from a system
    described by one of the following three forms:

        [1] x' = F(x, t, r, p)

        [2] M(x, p) x' = F(x, t, r, p)

        [3] M(q, p) u' = F(q, u, t, r, p)
            q' = G(q, u, t, r, p)

    where

        x : states
        t : time
        r : specified inputs
        p : constants
        q : generalized coordinates
        u : generalized speeds
        M : mass matrix (full or minimum)
        F : right hand side (full or minimum)
        G : right hand side of the kinematical differential equations

    The generated function is of the form F(x, t, p) or F(x, t, r, p)
    depending on whether the system has specified inputs or not.

    Parameters
    ==========
    right_hand_side : SymPy Matrix, shape(n, 1)
        A column vector containing the symbolic expressions for the
        right hand side of the ordinary differential equations. If the
        right hand side has been solved for symbolically then x' = f(x,
        t, r, p), if it hasn't then you must supply the accompanying
        mass matrix M(x) * x' = f(x, t, r, p).
    coordinates : sequence of SymPy Functions
        The generalized coordinates.
    speeds : sequence of SymPy Functions
        The generalized speeds.
    constants : sequence of SymPy Symbols
        All of the constants present in the right hand side and mass
        matrix.
    mass_matrix : sympy.Matrix, shape(n, n), optional
        This can be either the "full" mass matrix in M(x) * x' = f(x, t,
        r, p) or the minimal mass matrix in M(q) * u' = f(u, q, t, r,
        p).
    coordinate_derivatives : sympy.Matrix, shape(m, 1), optional
        If the "minimal" mass matrix is supplied, then this matrix
        represents the right hand side of q' = g(q, t, p).
    specifieds : sequence of SymPy Functions
        The specified exogenous inputs to the ODEs. These should be
        functions of time.
    linear_sys_solver : string or function
        Specify either `numpy` or `scipy` to use the linear solvers
        provided in each package or supply a function that solves a
        linear system Ax=b with the call signature x = solve(A, b). If
        you need to use custom kwargs for the scipy solver, pass in a
        lambda function that wraps the solver and sets them.
    constants_arg_type : string
        The generated function accepts two different types of arguments for
        the numerical values of the constants: an ndarray of the constants
        values in the correct order or a dictionary mapping the constants
        sysmbols to the numerical values. If None, this is determined inside
        of the generated function and can cause a significant slow down for
        performance critical code. If you know apriori what arg types you
        need to support choose either ``array`` or ``dictionary``. Note that
        ``array`` is faster than ``dictionary``.
    specifieds_arg_type : string
        The generated function accepts three different types of arguments
        for the numerical values of the specifieds: an ndarray of the
        specifieds values in the correct order, a function that generates
        the correct ndarray, a dictionary mapping the specifieds symbols or
        tuples of thereof to floats, ndarrays, or functions. If None, this
        is determined inside of the generated function and can cause a
        significant slow down for performance critical code. If you know
        apriori what arg types you want to support choose either ``array``,
        ``function``, or ``dictionary``. The speed of each, from fast to
        slow, are ``array``, ``function``, ``dictionary``.
    generator : string or and ODEFunctionGenerator, optional
        The method used for generating the numeric right hand
        side. The string options are {'lambdify'|'theano'|'cython'} with
        'lambdify' being the default. You can also pass in a custom subclass
        of ODEFunctionGenerator.

    Returns
    =======
    rhs_func : function
        A function which evaluates the derivaties of the states. See the
        function's docstring for more details after generation.

    Notes
    =====
    The generated function still supports the pre-0.3.0 extra argument
    style, i.e. args = {'constants': ..., 'specified': ...}, but only if
    ``constants_arg_type`` and ``specifieds_arg_type`` are both set to None.
    This functionality is deprecated and will be removed in 0.4.0, so it's
    best to adjust your code to support the new argument types.  See the
    docstring for the generated function for more info.

    """
    generators = {'lambdify': LambdifyODEFunctionGenerator,
                  'cython': CythonODEFunctionGenerator,
                  'theano': TheanoODEFunctionGenerator}

    try:
        # See if user passed in a custom class.
        g = generator(right_hand_side, coordinates, speeds, constants,
                      mass_matrix=mass_matrix,
                      coordinate_derivatives=coordinate_derivatives,
                      specifieds=specifieds,
                      linear_sys_solver=linear_sys_solver,
                      constants_arg_type=constants_arg_type,
                      specifieds_arg_type=specifieds_arg_type)
    except TypeError:
        # See if user passed in a string.
        try:
            g = generators[generator](right_hand_side, coordinates, speeds,
                                      constants, mass_matrix=mass_matrix,
                                      coordinate_derivatives=coordinate_derivatives,
                                      specifieds=specifieds,
                                      linear_sys_solver=linear_sys_solver,
                                      constants_arg_type=constants_arg_type,
                                      specifieds_arg_type=specifieds_arg_type)
        except KeyError:
            msg = '{} is not a valid generator.'.format(generator)
            raise NotImplementedError(msg)
        else:
            return g.generate()
    else:
        return g.generate()
