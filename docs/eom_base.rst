==============================
Two Mass Spring Damper Example
==============================

This example walks through the creation of a two mass spring damper system,
manually entered to the `SymbolicSystem` class and numerically simulated using
`pydy.system.System`. ::

    >>> import matplotlib.pyplot as plt
    >>> from numpy import linspace
    >>> from pydy.system import System
    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> import sympy.physics.mechanics.system as system

The first step is to define the dynamic and constant symbols that will be used
in the system. ::

    >>> x1, x2, u1, u2 = dynamicsymbols('x1 x2 u1 u2')
    >>> m1, c1, k1 = symbols('m1 c1 k1')
    >>> m2, c2, k2 = symbols('m2 c2 k2')

Next step is to define the equations of motion in multiple forms:

[1] x' = F(x, t, r, p)

[2] M(x, p) x' = F(x, t, r, p)

[3] M(q, p) u' = F(q, u, t, r, p)
    q' = G(q, u, t, r, p) ::

    >>> dyn_implicit_mat = Matrix([[m1,  0],
    ...                            [0,  m2]])
    >>> dyn_implicit_rhs = Matrix([-(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
    ...                            k2*x1 - k2*x2 + c2*u1 - c2*u2])
    >>> comb_implicit_mat = Matrix([[1, 0,  0,  0],
    ...                             [0, 1,  0,  0],
    ...                             [0, 0, m1,  0],
    ...                             [0, 0,  0, m2]])
    >>> comb_implicit_rhs = Matrix([u1,
    ...                             u2,
    ...                             -(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
    ...                             k2*x1 - k2*x2 + c2*u1 - c2*u2])
    >>> kin_explicit_rhs = Matrix([u1, u2])
    >>> comb_explicit_rhs = comb_implicit_mat.LUsolve(comb_implicit_rhs)

Now the various forms of varible entry to `SymbolicSystem` will be initiated ::

    >>> coordinates = (x1, x2)
    >>> speeds = (u1, u2)
    >>> states = (x1, x2, u1, u2)

This example will also show `system`'s ability to calculate an arbitrary user
defined functions along with the states. Specifically the kinetic energy and
potential energy  of the whole system is going to be determined. ::

    >>> KE = 0.5 * (m1*u1**2 + m2*u2**2)
    >>> PE = symbols("PE")
    >>> out_eqns = {"kinetic energy": KE, PE: 0.5*k1*x1**2 + 0.5*k2*(x2-x1)**2}

The `SymbolicSystem` instances are now ready to be initialized. To showcase
possible input combinations, all three of the equations of motion forms will be
used. ::

    >>> symsystem1 = system.SymbolicSystem(coordinates, comb_explicit_rhs, 
    ...                                    speeds=speeds, output_eqns=out_eqns)
    >>> symsystem2 = system.SymbolicSystem(states, comb_implicit_rhs, 
    ...                                    mass_matrix=comb_implicit_mat,
    ...                                    num_coordinates=2)
    >>> symsystem3 = system.SymbolicSystem(states, dyn_implicit_rhs, 
    ...                                    mass_matrix=dyn_implicit_mat,
    ...                                    coordinate_derivatives=kin_explicit_rhs)

There are some potential problems that can occur during the initialization of
the system.SymbolicSystem class. First, the class assumes if it recieves
coordinate_derivatives that the input equations of motion are in form [3]. This
means it will produce an error if a mass matrix is not recieved. Similarly if a
mass matrix is recieved and coordinate_derivatives are not specified, the class
is assuming that the equations of motion have been put in in form [2]. If the
equations of motion are input in form [1], the class will not have a notion of
the implicit form of the equatioins of motion and as such will return errors if
there is an attempt to access these attributes. ::

    >>> symsystem4 = system.SymbolicSystem(states, f, coordinate_derivatives=G)
    SyntaxError: Need to specify a mass matrix for eom form [3]
    >>> symsystem1.dyn_implicit_mat
    AttributeError: dyn_implicit_mat is not specified for eom form [1]
    >>> symsystem1.comb_implicit_mat
    AttributeError: comb_implicit_mat is not specified for eom form [1]
    >>> symsystem1.dyn_implicit_rhs
    AttributeError: dyn_implicit_rhs is not specified for eom form [1]
    >>> symsystem1.comb_implicit_rhs
    AttributeError: comb_implicit_rhs is not specified for eom form [1]

The system instance will now be initialized and set up to perform simulation. ::

    >>> sys = System(symsystem1)
    >>> sys.times = linspace(0, 10, num=100)
    >>> sys.constants = {m1: 10.0, c1: 10.0, k1: 10.0, m2: 5.0, c2: 5.0, k2:
    ...                  5.0}
    >>> sys.initial_conditions = {x1: 1.0, x2: 0.0, u1: 0.0, u2: 0.0}

Now the actual system simulation is ready to be run. ::

    >>> out = sys.integrate()
    >>> out
    array([[1.0, 0.0, 0.0, 0.0], [...], ...])

The user specified output equations should be able to be determined now that the
time simulation of the states is complete. These equations are calculated by
sys.output_eqns and so if there is an attempt to see the results first an error
will be returned. ::

    >>> symsystem1.output_eqns
    {"kinetic energy": KE, PE: 0.5*k1*x1**2 + 0.5*k2*(x2-x1)**2}
    >>> sys.output_eqns_calculate()
    {"kinetic energy": array([[0.0], [...], ...]),
     PE: array([[5.0], [...], ...])}
    >>> sym.output_eqns_results
    {"kinetic energy": array([[0.0], [...], ...]),
     PE: array([[5.0], [...], ...])}

With the simulation completed the output trajectories can be plotted using
matplotlib. ::

    >>> plt.plot(sys.times, out[:, 1])  
    >>> plt.plot(sys.times, out[:, 2])
    >>> plt.show()

System also has multiple plotting capabilities built into the class. ::

    >>> sys.plot_coordinates()
    >>> sys.plot_speeds()
    >>> sys.plot_states()
    >>> sys.plot_trajectories(x1, x2, u2, "kinetic energy", c1*x1+c2*x2)

The method `plot_trajectories()` can take as input different symbols contained
in the system, expressions using symbols defined in the system or keys of the
output equations dictionary.

=========================================
Simple Pendulum (x,y) Coordinates Example
=========================================

This code will go over the manual input of the equations of motion for the
simple pendulum into `Symbolic System` using x and y coordinates instead of
theta.

The equations of motion are formed at
http://nbviewer.jupyter.org/github/bmcage/odes/blob/master/docs/ipython/Planar%20Pendulum%20as%20DAE.ipynb` ::

    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> import sympy.physics.mechanics.system as system

The first step will be to initialize all of the dynamic and constant symbols. ::

    >>> x, y, u, v, lam = dynamicsymbols('x y u v lambda')
    >>> m, l, g = symbols('m l g')

Next step is to define the equations of motion in multiple forms:

[1] x' = F(x, t, r, p)

[2] M(x, p) x' = F(x, t, r, p)

[3] M(q, p) u' = F(q, u, t, r, p)
    q' = G(q, u, t, r, p) ::

    >>> dyn_implicit_mat = Matrix([[1, 0, -x/m],
    ...                            [0, 1, -y/m],
    ...                            [0, 0, l**2/m]])
    >>> dyn_implicit_rhs = Matrix([0, 0, u**2 + v**2 - g*y])
    >>> comb_implicit_rhs = Matrix([[1, 0, 0, 0, 0],
    ...                             [0, 1, 0, 0, 0],
    ...                             [0, 0, 1, 0, -x/m],
    ...                             [0, 0, 0, 1, -y/m],
    ...                             [0, 0, 0, 0, l**2/m]])
    >>> comb_implicit_rhs = Matrix([u, v, 0, 0, u**2 + v**2 - g*y])
    >>> kin_explicit_rhs = Matrix([u, v])
    >>> comb_explicit_rhs = comb_implicit_mat.LUsolve(comb_implicit_rhs)

Now the reference frames, points and particles will be set up so this
information can be passed into `system.SymbolicSystem` in the form of a bodies
and loads iterable. ::

    >>> N = ReferenceFrame('N')
    >>> A = N.orientnew('A', 'Axis', [theta, N.z])
    >>> O = Point('O')
    >>> P = O.locatenew('P', l * A.x)
    >>> Pa = Particle('Pa', P, m)

Now the bodies and loads iterables need to be initialized. ::

    >>> bodies = [Pa]
    >>> loads = [(P, g * m * N.x)]

The equations of motion are in the form of a differential algebraic equation
(DAE) and DAE solvers need to know which of the equations are the algebraic
expressions. This information is passed into `SymbolicSystem` as a list
specifying which rows are the algebraic equations. In this example it is a
different row based on the chosen equations of motion format. ::

    >>> alg_con = [2]
    >>> alg_con_full = [4]

An iterable containing the states now needs to be created for the solvers. ::

    >>> states = (x, y, u, v, lam)

Now the equations of motion instances can be created using the above mentioned
equations of motion formats. ::

    >>> symsystem1 = system.SymbolicSystem(states, comb_explicit_rhs, 
    ...                                    alg_con=alg_con_full, bodies=bodies, 
    ...                                    loads=loads)
    >>> symsystem2 = system.SymbolicSystem(states, comb_implicit_rhs, 
    ...                                    mass_matrix=comb_implicit_mat,
    ...                                    alg_con=alg_con_full,
    ...                                    coord_idxs=(0, 1))
    >>> symsystem3 = system.SymbolicSystem(states, dyn_implicit_rhs, 
    ...                                    mass_matrix=dyn_implicit_mat,
    ...                                    coordinate_derivatives=kin_explicit_rhs,
    ...                                    alg_con=alg_con, coord_idxs=(0, 1), 
    ...                                    speed_idxs=(2, 3))

The `SymbolicSystem` class can determine which of the states are considered
coordinates or speeds by passing in the indexes of the coordinates and speeds.
If these indexes are not passed in the object will not be able to differentiate
between coordinates and speeds. ::

    >>> symsystem1.coordinates
    AttributeError: The coordinates were not specified
    >>> symsystem2.speeds
    AttributeError: The speeds were not specified

Lastly here are some attributes accessible from the `SymbolicSystem` class. ::

    >>> symsystem1.states
    (x, y, u, v, lam)
    >>> symsystem2.coordinates
    (x, y)
    >>> symsystem3.speeds
    (u, v)
    >>> symsystem1.comb_explicit_rhs
    Matrix([[u(t)], [v(t)], [(-g*y(t) + u(t)**2 + v(t)**2)*x(t)/l**2],
            [(-g*y(t) + u(t)**2 + v(t)**2)*y(t)/l**2], [m*(-g*y(t) + u(t)**2 +
             v(t)**2)/l**2]])
    >>> symsystem2.comb_implicit_rhs
    Matrix([u, v, 0, 0, u**2 + v**2 - g*y])
    >>> symsystem2.comb_implicit_mat
    Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, -x/m], [0, 0, 0, 1,
             -y/m], [0, 0, 0, 0, l**2/m]])
    >>> symsystem3.dyn_implicit_rhs
    Matrix([0, 0, u**2 + v**2 - g*y])
    >>> symsystem3.dyn_implicit_mat
    Matrix([[1, 0, -x/m], [0, 1, -y/m], [0, 0, l**2/m]])
    >>> symsystem3.kin_explicit_rhs
    Matrix([u, v])
    >>> symsystem1.alg_con
    [4]
    >>> symsystem1.dynamic_symbols
    (x, y, u, v, lam)
    >>> symsystem1.constant_symbols
    (m, l, g)

Like coordinates and speeds, the bodies and loads attributes can only be
accessed if they are specified during initialization of the `SymbolicSystem`
class. ::

    >>> symsystem2.bodies
    AttributeError: The bodies were not specified
    >>> symsystem2.loads
    AttributeError: The loads were not specified

Several of the attributes are properties and as such do not support assignment.
These attributes are given below. ::

    >>> symsystem3.bodies = 42
    TypeError: Bodies does not support assignment
    >>> symsystem3.coordinates = 42
    TypeError: Coordinates does not support assignment
    >>> symsystem3.dyn_implicit_rhs = 42
    TypeError: dyn_implicit_rhs does not support assignment
    >>> symsystem3.comb_implicit_rhs
    TypeError: comb_implicit_rhs does not support assignment
    >>> symsystem3.loads = 42
    TypeError: Loads does not support assignment
    >>> symsystem3.dyn_implicit_mat = 42
    TypeError: dyn_implicit_mat does not support assignment
    >>> symsystem3.comb_implicit_mat = 42
    TypeError: comb_implicit_mat does not support assignment
    >>> symsystem3.kin_explicit_rhs = 42
    TypeError: kin_explicit_rhs does not support assignment
    >>> symsystem3.comb_explicit_rhs = 42
    TypeError: comb_explicit_rhs does not support assignment
    >>> symsystem3.speeds = 42
    TypeError: Speeds does not support assignment
    >>> symsystem3.states = 42
    TypeError: States does not support assignment
    >>> symsystem3.dynamic_symbols = 42
    TypeError: dynamic_symbols does not support assignment
    >>> symsystem1.constant_symbols = 42
    TypeError: constant_symbols does not support assignment


========================================
Simple Pendulum Theta Coordinate Example
========================================

This example walks through the same dynamical setup as ther previous but
defines the system by the angle theta instead of using x and y coordinates.
This results in an ODE system for the equations of motion rather than a DAE
system. Also the equations of motion will be formed by `LagrangesMethod` class
rather than being input manually. ::

    >>> from sympy import *
    >>> from sympy.physics.mechanics import LagrangesMethod, Lagrangian
    >>> from sympy.physics.mechanics import ReferenceFrame, Particle, Point
    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> from pydy.system import System

The first step is to create the dynamic and constant symbols used by the
system. ::

    >>> theta = dynamicsymbols('theta')
    >>> thetad = dynamicsymbols('theta', 1)
    >>> m, l, g = symbols('m l g')

Now the reference frames need to be set up. Reference frame A is set in the
plane perpendicular to the page containing segment OP. ::

    >>> N = ReferenceFrame('N')
    >>> A = N.orientnew('A', 'Axis', [theta, N.z])

The next step is to initialize the points and particles that will be used in
the dynamical system. ::

    >>> O = Point('O')
    >>> P = O.locatenew('P', l * A.x)
    >>> Pa = Particle('Pa', P, m)

With the points and reference frames determined, it is time to define how they
all move with respect to one another. ::

    >>> A.set_ang_vel(N, thetad * N.z)
    >>> O.set_vel(N, 0)
    >>> P.v2pt_theory(O, N, A)

Now the lagrangian and force list can be created and with these an instance of
`LagrangesMethod` can be initialized. ::

    >>> L = Lagrangian(N, Pa)
    >>> fl = [(P, g * m * N.x)]
    >>> l = LagrangesMethod(L, [theta], forcelist=fl, frame=N)

The `LagrangesMethod` instance can pass an instance of `SymbolicSystem` using
its `.to_system()` method. This allows the class to handle all of the
formatting for `SymbolicSystem` rather than making the user pass everything in
manually. For instance it will automatically change the equations to first
order form. ::

    >>> symsystem = l.to_system()
    >>> sys = System(symsystem)

Now that the system is set up, a simple time simulation will be performed. ::

    >>> sys.times = linspace(0, 10, num=100)
    >>> sys.constants = {m: 10, l: 5, g: 9.8}
    >>> sys.initial_conditions = {theta: 60, thetad: 0}
    >>> sys.integrate()
    array([[60.0, 0.0], [...], ...])

Display the kinetic energy change in time (obtained from the particle in the
bodies list). The kinetic energies are displayed in the order listed in the
`bodies` list. The last column is the kinetic energy of the whole system and is
simply the addition of all the other kinetic energies in the array at each time
step. ::

    >>> sys.body_kinetic_energies()
    array([[0.0, 0.0], [...], ...])

Here are some additional attributes accessible from the `SymbolicSystem` 
class. ::

    >>> symsystem.bodies
    [Pa]
    >>> symsystem.loads
    [(P, g * m * N.x)]
