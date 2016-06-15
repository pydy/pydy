==============================
Two Mass Spring Damper Example
==============================

This example walks through the creation of a two mass spring damper system,
manually entered to the `eombase` class and numerically simulated using
`pydy.system.System`. ::

    >>> import matplotlib.pyplot as plt
    >>> from numpy import linspace
    >>> from pydy.system import System
    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.mechanics import dynamicsymbols, eombase

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

    >>> mm = Matrix([[m1,  0],
    ...             [0,  m2]])
    >>> f = Matrix([-(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
    ...             k2*x1 - k2*x2 + c2*u1 - c2*u2])
    >>> mm_full = Matrix([[1, 0,  0,  0],
    ...                   [0, 1,  0,  0],
    ...                   [0, 0, m1,  0],
    ...                   [0, 0,  0, m2]])
    >>> f_full = Matrix([u1,
    ...                  u2,
    ...                  -(k1+k2)*x1 + k2*x2 - (c1+c2)*u1 + c2*u2,
    ...                  k2*x1 - k2*x2 + c2*u1 - c2*u2])
    >>> G = Matrix([u1, u2])
    >>> RHS = mm_full.LUsolve(f_full)

Now the various forms of varible entry to `eombase` will be initiated ::

    >>> coordinates = (x1, x2)
    >>> speeds = (u1, u2)
    >>> states = (x1, x2, u1, u2)

This example will also show `system`'s ability to calculate an arbitrary user
defined functions along with the states. Specifically the kinetic energy and
potential energy  of the whole system is going to be determined. ::

    >>> KE = 0.5 * (m1*u1**2 + m2*u2**2)
    >>> PE = symbols("PE")
    >>> out_eqns = {"kinetic energy": KE, PE: 0.5*k1*x1**2 + 0.5*k2*(x2-x1)**2}

The `eombase` instances are now ready to be initialized. To showcase possible
input combinations, all three of the equations of motion forms will be used. ::

    >>> eom1 = eombase.EOM(coordinates, RHS, speeds=speeds, 
    ...                    output_eqns=out_eqns)
    >>> eom2 = eombase.EOM(states, f_full, mass_matrix=mm_full,
    ...                    num_coordinates=2)
    >>> eom3 = eombase.EOM(states, f, mass_matrix=mm, coordinate_derivatives=G)

The system instance will now be initialized and set up to perform simulation. ::

    >>> sys = System(eom1)
    >>> sys.times = linspace(0, 10, num=100)
    >>> sys.constants = {m1: 10.0, c1: 10.0, k1: 10.0, m2: 5.0, c2: 5.0, k2:
    ...                  5.0}
    >>> sys.initial_conditions = {x1: 1.0, x2: 0.0, u1: 0.0, u2: 0.0}

Now the actual system simulation is ready to be run. ::

    >>> out = sys.integrate()
    >>> out
    array([[1.0, 0.0, 0.0, 0.0], [...], ...])

The user specified output equations should be able to be determine now that the
time simulation of the states is complete. ::

    >>> sys.output_eqns()
    {"kinetic energy": array([[0.0], [...], ...]),
     PE: array([[5.0], [...], ...])}
    >>> eom1.output_eqns
    {"kinetic energy": KE, PE: 0.5*k1*x1**2 + 0.5*k2*(x2-x1)**2}
    >>> eom1.output_eqns_results
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
simple pendulum into eombase using x and y coordinates instead of theta.

The equations of motion are formed at
http://nbviewer.jupyter.org/github/bmcage/odes/blob/master/docs/ipython/Planar%20Pendulum%20as%20DAE.ipynb` ::

    >>> from sympy import symbols, Matrix
    >>> from sympy.physics.mechanics import dynamicsymbols, eombase

The first step will be to initialize all of the dynamic and constant symbols. ::

    >>> x, y, u, v, lam = dynamicsymbols('x y u v lambda')
    >>> m, l, g = symbols('m l g')

Next step is to define the equations of motion in multiple forms:

[1] x' = F(x, t, r, p)

[2] M(x, p) x' = F(x, t, r, p)

[3] M(q, p) u' = F(q, u, t, r, p)
    q' = G(q, u, t, r, p) ::

    >>> mm = Matrix([[1, 0, -x/m],
    ...              [0, 1, -y/m],
    ...              [0, 0, l**2/m]])
    >>> f = Matrix([0, 0, u**2 + v**2 - g*y])
    >>> mm_full = Matrix([[1, 0, 0, 0, 0],
    ...                   [0, 1, 0, 0, 0],
    ...                   [0, 0, 1, 0, -x/m],
    ...                   [0, 0, 0, 1, -y/m],
    ...                   [0, 0, 0, 0, l**2/m]])
    >>> f_full = Matrix([u, v, 0, 0, u**2 + v**2 - g*y])
    >>> G = Matrix([u, v])
    >>> RHS = mm_full.LUsolve(f_full)

`Define bodies and loads`
Now the reference frames, points and particles will be set up so this
information can be passed into `eombase.EOM` in the form of a bodies and loads
iterable. ::

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
expressions. This information is passed into `eombase` as a list specifying
which rows are the algebraic equations. In this example it is a different row
based on the chosen equations of motion format. ::

    >>> alg_con = [2]
    >>> alg_con_full = [4]

An iterable containing the states now needs to be created for the solvers. ::

    >>> states = (x, y, u, v, lam)

Now the equations of motion instances can be created using the above mentioned
equations of motion formats. ::

    >>> eom1 = eombase.EOM(states, RHS, alg_con=alg_con_full, bodies=bodies,
    ...                    loads=loads)
    >>> eom2 = eombase.EOM(states, f_full, mass_matrix=mm_full,
    ...                    alg_con=alg_con_full, num_coordinates=2)
    >>> eom3 = eombase.EOM(states, f, mass_matrix=mm, coordinate_derivatives=G,
    ...                    alg_con=alg_con, num_coordinates=2, num_speeds=2)

Lastly here are some attributes accessible from the `EOM` class. ::

    >>> eom1.states
    (x, y, u, v, lam)
    >>> eom2.coordinates
    (x, y)
    >>> eom3.speeds
    (u, v)
    >>> eom1.rhs
    Matrix([[u(t)], [v(t)], [(-g*y(t) + u(t)**2 + v(t)**2)*x(t)/l**2], [(-g*y(t) + u(t)**2 + v(t)**2)*y(t)/l**2], [m*(-g*y(t) + u(t)**2 + v(t)**2)/l**2]])
    >>> eom2.forcing_full
    Matrix([u, v, 0, 0, u**2 + v**2 - g*y])
    >>> eom2.mass_matrix_full
    Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, -x/m], [0, 0, 0, 1, -y/m], [0, 0, 0, 0, l**2/m]])
    >>> eom3.forcing
    Matrix([0, 0, u**2 + v**2 - g*y])
    >>> eom3.mass_matrix
    Matrix([[1, 0, -x/m], [0, 1, -y/m], [0, 0, l**2/m]])
    >>> eom1.dynamic_symbols()
    (x, y, u, v, lam)
    >>> eom1.constant_symbols()
    (m, l, g)

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

The `LagrangesMethod` instance can pass an instance of eombase using its
`.to_eom()` method. This allows the class to handle all of the formatting for
eombase rather than making the user pass everything in manually. For instance
it will automatically change the equations to first order form. ::

    >>> eom = l.to_eom()
    >>> sys = System(EOM)

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

Here are some additional attributes accessible from the `eombase.EOM` class. ::

    >>> eom.bodies
    [Pa]
    >>> eom.loads
    [(P, g * m * N.x)]
