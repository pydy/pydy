.. jupyter-execute::

   #====================
   # rolling_ball_uneven_street
   #====================

.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`rolling_ball_uneven_street` or Jupyter notebook:
   :jupyter-download:notebook:`rolling_ball_uneven_street`.


.. jupyter-execute::

    import sympy as sm
    import sympy.physics.mechanics as me
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve, minimize
    import matplotlib.pyplot as plt
    import matplotlib as mp
    from mpl_toolkits.mplot3d import axes3d
    %matplotlib inline
    from IPython.display import HTML
    mp.rcParams['animation.embed_limit'] = 2**128 # I need this for large animations on my iPad
    from matplotlib import animation
    
    import pythreejs as p3js

A homogeneous, solid ball with radius r and mass m is running down an
uneven street without slipping. The ball is not allowed to jump, but is
in contact with the street at all times. The reaction forces to hold it
on the street are calculated.

The basic shape of the street is modeled by strassen_form, the
superimposed unevenness is modeled by strasse. (Strasse = Street in
German, gesamt = combined in German), that is gesamt = strasse +
strassen_form.

An ‘observer’, a particle of mass mo, may be attached anywhere inside
the ball.

For reasons which are not really clear to me, the total energy is not
constant, unless the ‘schmiegekreis’ (osculating circle as per google
translate) >> radius of the ball. There may be stiffness issues, as
method = ‘Radau’ in solve_ivp gives more reasonable results than no
method.

=================================================================================================================================
Define the sympy symbols and the dynamic symbols needed for the
calculation - q1, q2, q3 are the angle of the ball w.r.t. N. - u1, u2,
u3 are the speeds - x, z are the coordinates of the street - ux, uz the
speeds. - auxx, auxy, auxz, fx, fy, fz are needed for the reaction
forces at the contact point CP - the rest are constants needed to
describe the physical properties of the system: - reibung is the
friction - amplitude and frequenz determine the shape of the street - $
:raw-latex:`\alpha`, :raw-latex:`\beta`, :raw-latex:`\gamma `$ determine
the location of the observer w.r.t. the center of the ball.

.. jupyter-execute::

    start = time.time()
    
    q1, q2, q3 = me.dynamicsymbols('q1 q2 q3') 
    u1, u2, u3 = me.dynamicsymbols('u1 u2 u3')
    
    x, z = me.dynamicsymbols('x z')            
    ux, uz = me.dynamicsymbols('ux uz')
    
    auxx, auxy, auxz, fx, fy, fz = me.dynamicsymbols('auxx, auxy, auxz, fx, fy, fz')
    
    
    m, mo, g, r, iXX, iYY, iZZ = sm.symbols('m, mo, g, r, iXX, iYY, iZZ')
    amplitude, frequenz, reibung, alpha, beta, gamma, t = sm.symbols('amplitude frequenz reibung alpha beta gamma t')


Define the frames and the points as follows: - N = inertial frame - A2 =
frame attached to the ball - A3 = auxiliary frame to align the frame
used here with the frame (I think) needed for Pythreejs - P0 = fixed
reference point - CP = contact point, where the ball touches the street.
- Dmc = center of the ball - m_Dmc = particle anywhwere within the ball

.. jupyter-execute::

    N, A2, A3 = sm.symbols('N, A2, A3', cls = me.ReferenceFrame)         
    P0, CP, Dmc, m_Dmc = sm.symbols('P0, CP, Dmc, m_Dmc', cls = me.Point)
    
    P0.set_vel(N, 0.)


Here the street is modelled. The larger rumpel, the more ‘uneven’ the
street will be. The smaller frequenz, the ‘flatter’ the street. However
larger rumple values will increase the force terms. rumpel must be an
integer.

.. jupyter-execute::

    #============================================
    rumpel = 2                          
    #============================================
    strasse = sum([amplitude/j * (sm.sin(j*frequenz * x) + sm.sin(j*frequenz * z)) for j in range(1, rumpel)])
    strassen_form = (frequenz/2. * x)**2  + (frequenz/2. * z)**2
    gesamt = strassen_form  + strasse


If Jakob is set to True, the jacobian will be calculated and used in
solve_ivp. In my test runs, it did not help any. Maybe this is so
because it has over 900,000 operations and the numerical errors ‘kill’
any effect the jacobian may have.

.. jupyter-execute::

    #============================================
    Jakob = False                       
    #============================================


This is to determine the maximum radius of the ball so anywhere on the
street it has only one contact point. the radius of the ball must be
smaller than the smallest osculating circle (Schmiegekreis) of the road.
Here I calculate this for the X and for the Z direction and simply take
the smaller one of the two. (I did not bother to try and find the
correct formula for higher dimensional functions)

.. jupyter-execute::

    r_max_x = (sm.S(1.) + (gesamt.diff(x))**2 )**sm.S(3/2)/gesamt.diff(x, 2)
    r_max_z = (sm.S(1.) + (gesamt.diff(z))**2 )**sm.S(3/2)/gesamt.diff(z, 2)

Get the frames oriented. $ rot $ and $ rot1 $ are used below for the
kinematic equations of motion. Frame A3 is needed to align the
orientation used to set up Kane’s equation, with the orientation (I
think) needed for pythreejs animation.

NOTE: in order to get the ‘smallest’ equations of motion, it is
necessary to use A2, not N, in defining the speed. That is A2.x, not
N.x, etc. in;

A2.set_ang_vel(N, u1 \* A2.x + u2 \* A2.y + u3 \* A2.z).

This is related to the ‘inner workings’ of orient_body_fixed. See also
the comment further down where the kinematic equations are formed.

.. jupyter-execute::

    A2.orient_body_fixed(N, (q1, q2, q3), '123') 
    rot = (A2.ang_vel_in(N))
    A2.set_ang_vel(N, u1*A2.x + u2*A2.y + u3*A2.z)
    rot1 = (A2.ang_vel_in(N))
    
    A3.orient_axis(A2, -sm.pi/2., A2.z)
    A3.set_ang_vel(A2, 0.)

CP is the contact point, where the ball touches the street. As it is
part of the street, it has no real velocity. The
:math:`auxx, auxy, auxz` are ‘virtual’ speeds needed for the reaction
forces :math:`fx, fy, fz`.

.. jupyter-execute::

    CP.set_pos(P0, x*N.x + gesamt*N.y + z*N.z)
    CP.set_vel(N, auxx*N.x + auxy*N.y + auxz*N.z) 

vektor is the normal to the streel at point (x, z), the X/Z coordinates
of the contact point CP. I found this formula in the internet. It points
‘inwards’, hence the leading minus sign. Dmc is in the direction of this
normal at distance r from CP

.. jupyter-execute::

    vektor = -(gesamt.diff(x)*N.x - N.y + gesamt.diff(z)*N.z)
    Dmc.set_pos(CP, r * vektor.normalize())

The instantaneous contact point CP is fixed in A2 - for that instant.
Hence, v2pt_theory is correct for DMC. (The semicolon at the very end is
only to avoid that the speed of m_Dmc will be printed)

.. jupyter-execute::

    Dmc.v2pt_theory(CP, N, A2) 
    
    m_Dmc.set_pos(Dmc, r*(alpha*A2.x + beta*A2.y + gamma*A2.z))
    m_Dmc.v2pt_theory(Dmc, N, A2);

The instantaneous contact point has no speed, as it is part of the
street, too. We need :math:`x(t), z(t)` to get the location of
subsequent contact points. These subsequent contact points will be
located at :math:`(x(t), gesamt(x(t), z(t)), z(t))` .

Relationship of x(t) to q(t): Obviously, $ x(t) = function(q(t),
gesamt(x(t), z(t)), r) $. When the ball is rotated through an angle
:math:`q`, the arc length is :math:`r * q(t)`.

The arc length of a function f(k(t)) from 0 to :math:`x(t)` is: $
:raw-latex:`\int`\_{0}^{x(t)} :raw-latex:`\sqrt`(1 + d/dx(f(k(t)))^2)
,dk  $ (I found this in the internet)

This gives the sought after relationship between :math:`q(t)` and
:math:`x(t)`: $ r \* (-q(t)) = :raw-latex:`\int`\_{0}^{x(t)}
:raw-latex:`\sqrt`(1 + d/dk(gesamt(k(t), z(t)))^2) ,dk  $,
differentiated: - $ r \* (-u) = :raw-latex:`\sqrt`(1 + d/dx(gesamt(x(t),
z(t)))^2) \* d/dt(x(t)) $, that is solved for d/dt(x(t)): - $ d/dt(x(t))
= :raw-latex:`\frac{-(r * u)}` {:raw-latex:`\sqrt`(1 + d/dx(gesamt(x(t),
z(t)))^2))} $

A similar formula holds for the speed in Z direction

$ d/dt(x(t)) = :raw-latex:`\frac{- (u3 * r)}` {:raw-latex:`\sqrt`(1 +
d/dx(gesamt(x(t), z(t)))^2)} $

$ d/dt(z(t)) = :raw-latex:`\frac{(u1 * r)}` {:raw-latex:`\sqrt`(1 +
d/dz(gesamt(x(t), z(t))))^2)} $

As speeds are vectors, they can be added to give the resultant speed.

The + / - signs are a consequence of the ‘right hand rule’ for frames.
These are the sought after first order differential equations for
:math:`(x(t), z(t))`.

.. jupyter-execute::

    OMEGA = u1*A2.x + u2*A2.y + u3*A2.z
    u3_wirk = me.dot(OMEGA, N.z)
    u1_wirk = me.dot(OMEGA, N.x)
    rhsx = (-u3_wirk * r / sm.sqrt(1. + gesamt.diff(x)**2))
    rhsz =  (u1_wirk * r / sm.sqrt(1. + gesamt.diff(z)**2))

Now, the standard machinery to set up Kane’s equation starts.

As the energy of a system may give hints that maybe a mistake was made
setting up Kane’s equations, e.g. absebt any friction or applied
forces/torques the total energy must be constant, I like to calculate
them and look at them.

.. jupyter-execute::

    I = me.inertia(A2, iXX, iYY, iZZ)                                              
    Body = me.RigidBody('Body', Dmc, A2, m, (I, Dmc))
    observer = me.Particle('observer', m_Dmc, mo)
    BODY = [Body, observer]
    
    subs_dict = {sm.Derivative(x, t): rhsx, sm.Derivative(z, t): rhsz}
    energie_dict = {i: 0. for i in (auxx, auxy, auxz, fx, fy, fz)}
    pot_energie = (m * g * me.dot(Dmc.pos_from(P0), N.y) + mo * g * me.dot(m_Dmc.pos_from(P0), N.y)
                  ).subs(subs_dict).subs(energie_dict)
    kin_energie = (Body.kinetic_energy(N) + observer.kinetic_energy(N)).subs(subs_dict).subs(energie_dict)

The forces :math:`fx, fy, fz` are the reaction forces on CP, the contact
point

.. jupyter-execute::

    FL1 = [(Dmc, -m*g*N.y), (m_Dmc, -mo*g*N.y), (CP, fx*N.x + fy*N.y + fz*N.z)]
    Torque = [(A2, -reibung * (u1*A2.x + u2*A2.y + u3*A2.z))]
    FL = FL1  + Torque

It is very important to use A2 when forming the kinematic equations of
motion. Using N will increase the number of operations in force and in
eingepraegt by a factor of around 7…8 in this case. I have seen worse
differences.

.. jupyter-execute::

    kd = [me.dot(rot - rot1, uv) for uv in A2]
    
    q = [q1, q2, q3]
    u = [u1, u2, u3]
    aux = [auxx, auxy, auxz]
    
    KM = me.KanesMethod(N, q_ind=q, u_ind=u, kd_eqs=kd, u_auxiliary=aux)
    (fr, frstar) = KM.kanes_equations(BODY, FL)
    MM = KM.mass_matrix_full
    force = KM.forcing_full

| Add :math:`rhsx, rhsz` at the bottom of force, to get
  :math:`d/dt(x) = rhsx`, same for z. This is to numerically integrate
  x(t), z(t). Remember:
| - $ d/dt(x(t)) = rhsx $
| - $ d/dt(z(t)) = rhsz $

I find it useful to print the dynamic symbols and the free symbols of
such functions. This may give hints that something was wrong in
describing the system. The number of operations is also useful to know,
I think.

.. jupyter-execute::

    force = (sm.Matrix.vstack(force, sm.Matrix([rhsx, rhsz]))).subs(subs_dict)
    print('force DS', me.find_dynamicsymbols(force))
    print('force free symbols', force.free_symbols)
    print('force has {} operations'.format(sum([force[i].count_ops(visual=False) 
                    for i in range(len(force))])), '\n')

The mass matrix MM must be enlarged accordingly

CP_pos, Dmc_pos, m_Dmc_pos, vektor_loc are needed later for plotting
only.

.. jupyter-execute::

    MM = sm.Matrix.hstack(MM, sm.zeros(6, 2))
    hilfs = sm.Matrix.hstack(sm.zeros(2, 6), sm.eye(2))
    MM = sm.Matrix.vstack(MM, hilfs )
    print('MM DS', me.find_dynamicsymbols(MM))
    print('MM free symbols', MM.free_symbols)
    print('MM has {} operations'.format(sum([MM[i, j].count_ops(visual=False) 
                    for i in range(MM.shape[0]) for j in range(MM.shape[1])])), '\n')
    
    
    CP_pos = [(me.dot(CP.pos_from(P0), uv)).subs(subs_dict) for uv in N]
    Dmc_pos = [(me.dot(Dmc.pos_from(P0), uv)).subs(subs_dict) for uv in N]# for later plotting only
    m_Dmc_pos = [(me.dot(m_Dmc.pos_from(P0), uv)).subs(subs_dict) for uv in N] #later plotting only
    
    vektor_loc = [me.dot(vektor, uv) for uv in N]

Find the reaction forces. As they depend on $ rhs = MM^{-1} \* force $
and rhs becomes quite large, it is calculated numerically later. Hence,
RHS ‘substitutes’ for rhs to be calculated numerically later

.. jupyter-execute::

    RHS = [sm.symbols('rhs' + str(i)) for i in range(8)]
    eingepraegt_dict = {sm.Derivative(i, t): RHS[j] for j, i in enumerate(q + u)}
    eingepraegt = ((KM.auxiliary_eqs).subs(eingepraegt_dict)).subs(subs_dict)
    print('eingepraegt DS', me.find_dynamicsymbols(eingepraegt))
    print('eingepraegt free symbols', eingepraegt.free_symbols)
    print('eingepraegt has {} operations'.format(sum([eingepraegt[i].count_ops(visual=False) 
                    for i in range(len(eingepraegt))])), '\n')

Form the Jacobian Matrix, this should help solve_ivp, but it does not
seem to do any good.

.. jupyter-execute::

    if Jakob == True:
        jakob = force.jacobian(sm.Matrix(q + u + [x, z]))
        print('jakob shape', jakob.shape)
        print('jakob DS', me.find_dynamicsymbols(jakob))
        print('jakob free symbols', jakob.free_symbols)
        print('jakob has {} operations'.format(sum([jakob[i, j].count_ops(visual=False) 
                    for i in range(jakob.shape[0]) for j in range(jakob.shape[1])])), '\n')


Convert all the sympy functions to numpy functions, so numerical
calculations may be done.

.. jupyter-execute::

    pL = [m, mo, g, r, iXX, iYY, iZZ, amplitude, frequenz, reibung] + [alpha, beta, gamma]
    qL = q + u + [x, z]
    F = [fx, fy, fz]
    
    MM_lam = sm.lambdify(qL + pL, MM, cse=True)
    force_lam = sm.lambdify(qL + pL, force, cse=True)
    
    CP_pos_lam = sm.lambdify(qL + pL, CP_pos, cse=True)
    Dmc_pos_lam = sm.lambdify(qL + pL, Dmc_pos, cse=True)
    m_Dmc_pos_lam = sm.lambdify(qL + pL, m_Dmc_pos, cse=True)
    
    gesamt_lam = sm.lambdify([x, z] + [amplitude, frequenz], gesamt, cse=True)
    
    pot_lam = sm.lambdify(qL + pL, pot_energie, cse=True)
    kin_lam = sm.lambdify(qL + pL, kin_energie, cse=True)
    
    eingepraegt_lam = sm.lambdify(F + qL + pL + RHS, eingepraegt, cse=True)
    if Jakob == True:
        jakob_lam = sm.lambdify(qL + pL, jakob, cse=True)
    
    r_max_lam = sm.lambdify([x, z] + pL, [r_max_x, r_max_z], cse=True)
    
    vektor_loc_lam = sm.lambdify(qL + pL, vektor_loc, cse=True)
    
    print('it took {:.3f} sec to establish Kanes equations'.format(time.time() - start))


Perform the numerical integration.

Input data: - m1: mass of the ball - mo1: mass of the particle
(observer) - amplitude1, frequenz1: parameters of the street (note: the
smaller frequenz1, the more even the street) -
:math:`\alpha, \beta, \gamma`: define the location of the particle
relative to the center of the ball - :math:`q11, q21, q31`: initial
generalized coordinates - :math:`u11, u21, u31`: initial speeds - x0,
z0: initial coordinates of the contact point - intervall: the
integration will run from 0 to intervall - schritte: number of time
instances returned by ivp_solve

.. jupyter-execute::

    #================================================================
    m1 = 1.e0
    mo1 = 1.e0
    r1 = 2.
    amplitude1 = 0.1
    frequenz1 = 1.e-1
    reibung1 = 0.
    
    alpha1, beta1, gamma1 = 0.95, 0., 0.
    
    q11, q21, q31 = 0., 0., 0.
    u11, u21, u31 = 5., 5., 5.
    x0, z0 = 0., 0.
    intervall = 25.
    schritte = 250
    
    #================================================================
    start = time.time()
    print('Arguments')
    print('[m, g, r, iXX, iYY, iZZ, amplitude, frequenz, reibung, alpha, beta, gamma]')

Ensure that the particle (observer) is inside the ball, that is $
:raw-latex:`\alpha`^2 + :raw-latex:`\beta`^2 + :raw-latex:`\gamma`^2 <
1. $ If this is not the case, an Exception is raised.

.. jupyter-execute::

    if alpha1**2 + beta1**2 + gamma1**2 >= 1.:
        raise Exception('center of mass outside of the ball')
    
    iXXe = 2./5. * m1 * r1**2
    
    pL_vals = [m1, mo1, 9.8, r1, iXXe, iXXe, iXXe, amplitude1, frequenz1, reibung1, alpha1, beta1, gamma1]
    print(pL_vals, '\n')
    
    y0 = [q11, q21, q31, u11, u21, u31] + [x0, z0]

The ball must touch the street at exactly one point only. The radius of
the ball r must be small than the smallest bevel radius. If the smallest
bevel curve radius (Schmiegekreis) is smaller than the radius of the
ball, an exception is raised.

.. jupyter-execute::

    def func1(x, args):
    # just needed to get the arguments matching for minimuze
        return np.abs(r_max_lam(*x, *args)[0])
    
    def func2(x, args):
    # just needed to get the arguments matching for minimuze
        return np.abs(r_max_lam(*x, *args)[1])
    
    x0 = (0.1, 0.1)      # initial guess
    minimal1 = minimize(func1, x0, pL_vals)
    minimal2 = minimize(func2, x0, pL_vals)
    
    minimal = min(minimal1.get('fun'), minimal2.get('fun'))
    
    if pL_vals[3] < minimal:
        print('selected radius = {} is less than minimally admissible radius = {:.2f}, hence o.k.'
              .format(pL_vals[3], minimal), '\n')
    else:
        print('selected radius {} is larger than admissible radius {:.2f}, hence NOT o.k.'
              .format(pL_vals[3], minimal), '\n')
        raise Exception('Radius of ball is too large')


Now the actual numerical integration starts.

.. jupyter-execute::

    times = np.linspace(0, intervall, schritte)
    
    if Jakob == True:
        def JAKOB(t, y, args):
            vals = np.concatenate((y, args))
            sol = np.linalg.solve(MM_lam(*vals), jakob_lam(*vals))
            return np.array(sol)
    else:
        JAKOB = None
        
    def gradient(t, y, args):
        vals = np.concatenate((y, args))
        sol = np.linalg.solve(MM_lam(*vals), force_lam(*vals))
        return np.array(sol).T[0]
    
    t_span = (0., intervall)
    resultat1 = solve_ivp(gradient, t_span, y0, t_eval=times, args=(pL_vals,), method='Radau'
                          , atol=1.e-6, rtol=1.e-3, jac=JAKOB)
    resultat = resultat1.y.T
    print(resultat.shape)
    event_dict = {-1: 'Integration failed', 0: 'Integration finished successfully', 1: 'some termination event'}
    print(event_dict[resultat1.status])
        
    print('resultat shape', resultat.shape, '\n')
    
    print("To numerically integrate an intervall of {} sec, with jacobi = {} the routine cycled {} times and it took {:.3f} sec "
          .format(intervall, Jakob, resultat1.nfev, time.time() - start))


Plot the rotational speeds of the ball.

.. jupyter-execute::

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, j in zip(range(3, 6), ('u1', 'u2', 'u3')):
        ax.plot(times, resultat[:, i], label = str(j))
    ax.set_title('rotational speed of ball')
    ax.legend();

Plot the energies of the ball. Total energy must be constant - which it
is not always, due to numerical issues in my opinion.

.. jupyter-execute::

    pot_np = np.empty(schritte)
    kin_np = np.empty(schritte)
    total_np = np.empty(schritte)
    for l in range(schritte):
        pot_np[l] = pot_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)
        kin_np[l] = kin_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)
        total_np[l] = pot_np[l] + kin_np[l]
    
    if reibung1 == 0.:
        fehler = ((x:=max(total_np)) - min(total_np)) / x * 100.
        print('max deviation from total energy = constant is {:.3f} % of max total energy'.format(fehler))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, kin_np, label = 'kin energy')
    ax.plot(times, pot_np, label = 'pot energy')
    ax.plot(times, total_np, label = 'total energy')
    ax.set_title('Energy of the system, friction = {}'.format(reibung1))
    ax.legend();


Plot the locations of CP and Dmc.

.. jupyter-execute::

    CPx = np.empty(schritte)
    CPy = np.empty(schritte)
    CPz = np.empty(schritte)
    
    gesamt_np = np.empty(schritte)
    
    Dmcx = np.empty(schritte)
    Dmcy = np.empty(schritte)
    Dmcz = np.empty(schritte)
    
    for l in range(schritte):
        CPx[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[0]
        CPy[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[1]
        CPz[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[2]
        
        Dmcx[l] = Dmc_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[0]
        Dmcy[l] = Dmc_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[1]
        Dmcz[l] = Dmc_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[2]
        
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, CPx, label = 'CPx')
    ax.plot(times, CPy, label = 'CPy')
    ax.plot(times, CPz, label = 'CPz')
    
    ax.plot(times, Dmcx, label = 'Dmcx')
    ax.plot(times, Dmcy, label = 'Dmcy')
    ax.plot(times, Dmcz, label = 'Dmcz')
    ax.set_title('Locations of various points')
    ax.legend();


Calculate and plot the reaction forces. 1. calculate $ rhs = MM^{-1} \*
force $ numerically, this is faster than doing it symbolically 2. find
the reaction forces numerically by solving :math:`eingepraegt = 0` for
:math:`fx, fy, fz` 3. Plot the results

.. jupyter-execute::

    RHS1 = np.zeros((schritte, resultat.shape[1]))
    for i in range(schritte):
        RHS1[i, :] = np.linalg.solve(MM_lam(*[resultat[i, j]for j in range(resultat.shape[1])], *pL_vals), 
            force_lam(*[resultat[i, j] for j in range(resultat.shape[1])], 
            *pL_vals)).reshape(resultat.shape[1])
    
    def func (x, *args):
    # just serves to 'modify' the arguments for fsolve.
        return eingepraegt_lam(*x, *args).reshape(3)
    
    kraftx = np.zeros(schritte)
    krafty = np.zeros(schritte)
    kraftz = np.zeros(schritte)
    x0 = tuple((1., 1., 1.))   # initial guess
    
    for i in range(schritte):
        y0 = [resultat[i, j] for j in range(resultat.shape[1])]
        rhs = [RHS1[i, j] for j in range(8)]
        args = tuple(y0 + pL_vals + rhs)
        A = fsolve(func, x0, args=args).reshape(3)
        x0 = tuple(A)      # improved guess. Should speed up convergence of fsolve
        kraftx[i] = A[0]
        krafty[i] = A[1]
        kraftz[i] = A[2]
    
    fig, ax = plt.subplots(figsize=(10, 5))
    plt.plot(times, kraftx, label='X force')
    plt.plot(times, krafty, label='Y force')
    plt.plot(times, kraftz, label='Z force')
    ax.set_title('Reaction Forces on CP')
    plt.legend();


This gives the path of the ball, projected on the X / Z horizontal
plane. The height above the lowest point of the street is indicated by
its color. This seems to give a better idea on how it rolls, compared to
the pythreejs animation below - at least to me.

.. jupyter-execute::

    CPx = np.empty(schritte)
    CPy = np.empty(schritte)
    CPy1 = np.empty(schritte)
    CPz = np.ones(schritte)
    for l in range(schritte):
        CPx[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[0]
        CPy[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[1]
        CPz[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[2]
    
    xmin1 = min(CPx)
    xmin2 = min(CPz)
    xmin = min(xmin1, xmin2)
    
    xmax1 = max(CPx)
    xmax2 = max(CPz)
    xmax = max(xmax1, xmax2)


This is to asign colors of ‘plasma’ to the points, and scale them
properly, so that the lowest point gets the darkest color of plasma, and
the highest point the ball reaches gets the brightest color.

.. jupyter-execute::

    for i in range(len(CPy)):
        CPy1[i] = int(CPy[i])
    ymin = min(CPy1)
    ymax = max(CPy1)
    
    Test = mp.colors.Normalize(ymin, ymax)
    Farbe = mp.cm.ScalarMappable(Test, cmap='plasma')
    farbe1 = Farbe.to_rgba(CPy1[0])    # color of the starting position

Now the animation itself starts. Note: the last command, HTML(…) is
needed to display the anmation on an iPad. It may not be needed for
another machine, I do not know

.. jupyter-execute::

    def animate_pendulum(times, x1, y1, z1):
        
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'aspect': 'equal'})
        fig.colorbar(Farbe, label='Height of the ball above ground level', shrink=0.9, ax=ax)
        ax.axis('on')
        ax.set_xlim(xmin - 1., xmax + 1.)
        ax.set_ylim(xmin - 1., xmax + 1.)
        ax.set_xlabel(' X Axis', fontsize=15)
        ax.set_ylabel('Z Axis', fontsize=15)
    
        line1, = ax.plot([], [], 'o', markersize=15, color=farbe1)       # starting point
        line2, = ax.plot([], [], 'o', markersize=15)                     # final point
        line3, = ax.plot([], [], color='blue', linewidth=0.25)           # to trace the movement of CP
        
        def animate(i):
            farbe2 = Farbe.to_rgba(y1[i])                                # color of the actual point at time i
            ax.set_title('running time {:.1f} sec'.format(times[i]), fontsize=15)
            line1.set_data(x1[0], z1[0])
            line2.set_data(x1[i], z1[i])
            line2.set_color(farbe2)
            line3.set_data(x1[: i], z1[: i])   
            return line1, line2, line3,
    
        anim = animation.FuncAnimation(fig, animate, frames=len(times),
                                       interval=1000*max(times) / len(times),
                                       blit=True)
        plt.close(fig)
        return anim
    
    anim = animate_pendulum(times, CPx, CPy, CPz)
    HTML(anim.to_jshtml())    


Animation using pythreejs.
This is basically copied from a program by Jason Moore, just adapted to my needs here.

NOTE: the 'reference frame' for pythreejs seems to be:
- X - axis downwards, color red
- Y - axis to the right, color green (hence:)
- Z - axis pointing to the observer, color blue

At the end of this program, sometimes some warning comes up. This seems to be harmless, everything works fine, at least
on my machine.

The 4 x 4 matirces TB, TC hold the components of the frame attached to the ball and the particle respectivel; in the upper left 3 x 3 sub matrix. (Here it is the same A2 frame, of course.)
The location of the origin of the frame is stored in the lower left 1 x 3 matrix. These locations of course are different for ball and particle.

.. jupyter-execute::

    xmax1 = max(CPx)
    xmax2 = max(CPz)
    xmax = max(xmax1, xmax2)
    
    
    TB = sm.eye(4)
    TB[:3, :3] = A3.dcm(N)
    TB = TB.reshape(16, 1)
    
    TC = sm.eye(4)
    TC[:3, :3] = A3.dcm(N)
    TC = TC.reshape(16, 1)
    
    TB_lam = sm.lambdify(qL + pL, TB) 
    TC_lam = sm.lambdify(qL + pL, TC)


Create the TBs / TCs, containing the information of the rotation for
every time step. This information is obviously a function of q1, q2, q3.

scala is an auxuíliary value to adjust the picture, I know no better
way.

TBi[12], TBi[13], TBi[14] hold the location Dmc, same for TCi and m_Dmc.
(the location of the origins of the respective body fixed frames, as
mentioned above)

As the axis chosen for solving the equations of motion, and the axis
given by pythreejs do not coincide, the values for TBi[..] / TCi[..]
must be permutated accordingly.

.. jupyter-execute::

    TBs = []   # for the ball
    TCs = []   # for its center ofgravity
    #=============================================================
    scala = 1.
    #=============================================================
    for k in range(resultat.shape[0]):
        TBi = TB_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)  # same for ball and for
        TCi = TC_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)  # the center of mass
        
        TBi[12] = -Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[1] /scala
        TBi[13] = Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[0] / scala
        TBi[14] = Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[2] / scala
        
        TCi[12] = -m_Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[1] / scala
        TCi[13] = m_Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[0] / scala
        TCi[14] = m_Dmc_pos_lam(*[resultat[k, l] for l in range(resultat.shape[1])], *pL_vals)[2] / scala
    
        TBs.append(TBi.squeeze().tolist())
        TCs.append(TCi.squeeze().tolist())

Create the objects, which will move 1. The ball 2. The particle

.. jupyter-execute::

    body_geom_B = p3js.SphereGeometry(r1, 12, 12)
    body_material_B = p3js.MeshStandardMaterial(color='orange', wireframe=True)
    body_mesh_B = p3js.Mesh(geometry=body_geom_B, material=body_material_B, name='ball')
    
    body_geom_C = p3js.SphereGeometry(r1/5., 12, 12)
    body_material_C = p3js.MeshStandardMaterial(color='blue', wireframe=False)
    body_mesh_C = p3js.Mesh(geometry=body_geom_C, material=body_material_C, name='zentrum')
    
    body_mesh_B.matrixAutoUpdate = False
    body_mesh_B.add(p3js.AxesHelper(r1*2))  # length of the axis of the ball system A2
    body_mesh_B.matrix = TBs[0]             # starting point of the animation
    
    body_mesh_C.matrixAutoUpdate = False
    body_mesh_C.add(p3js.AxesHelper(0.01))    # length of the axis of the center of mass system A2
    body_mesh_C.matrix = TCs[0]             # starting point of the animation


Create the ‘picture’.

All the ‘parameters’ are taken by trial and error, I am not very sure
what they mean.

.. jupyter-execute::

    view_width = 1200
    view_height = 400
    
    camera = p3js.PerspectiveCamera(position=[0.0, 0.0, 1.5 * xmax],
                                    up=[-1.0, 0.0, 0.0],
                                    aspect=view_width/view_height)
    
    key_light = p3js.DirectionalLight(position=[0, 0, 10])
    ambient_light = p3js.AmbientLight()
    
    axes = p3js.AxesHelper(20)
    
    children = [body_mesh_B, axes, camera, key_light, ambient_light] + [body_mesh_C, axes, camera, key_light, ambient_light]
    
    scene = p3js.Scene(children=children)
    
    controller = p3js.OrbitControls(controlling=camera)
    renderer = p3js.Renderer(camera=camera, scene=scene, controls=[controller],
                             width=view_width, height=view_height)
    
    
    track_B = p3js.VectorKeyframeTrack(
        name="scene/ball.matrix",
        times=times,
        values=TBs)
    
    track_C = p3js.VectorKeyframeTrack(
        name="scene/zentrum.matrix",
        times=times,
        values=TCs)
    
    
    
    tracks = [track_B, track_C]
    duration = times[-1] - times[0]
    clip = p3js.AnimationClip(tracks=tracks, duration=duration)
    action = p3js.AnimationAction(p3js.AnimationMixer(scene), clip, scene)
    renderer

Running this somehow starts everything. I do not know the details of
pythreejs

.. jupyter-execute::

    action

Plot approximate shape of the street. It does not really show the
unevenness, just the overall shape.

.. jupyter-execute::

    CPx = np.empty(schritte)
    CPy = np.empty(schritte)
    CPz = np.ones(schritte)
    for l in range(schritte):
        CPx[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[0]
        CPy[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[1]
        CPz[l] = CP_pos_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)[2]
        
    # needed to give the picture the right size.
    xmin1 = min(CPx)
    xmin2 = min(CPz)
    xmin3 = min(xmin1, xmin2)
    
    xmax1 = max(CPx)
    xmax2 = max(CPz)
    xmax3 = max(xmax1, xmax2)
    
    xmin = max(np.abs(xmin3), abs(xmax3))
    xs = np.linspace(-xmin, xmax, 500)
    ys = np.linspace(-xmin, xmax, 500)
    X, Y = np.meshgrid(xs, ys)
    Z = gesamt_lam(X, Y, amplitude1, frequenz1)
    
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    
    strasse = ax.plot_surface(X, Y, Z, cmap='plasma', linewidth=0., antialiased=True)
    ax.set_xlabel('X - Axis', fontsize = 15)
    ax.set_ylabel('Z - Axis', fontsize = 15)
    ax.set_title('Approximate shape of the street', fontsize=12)
    ax.set_zlabel('Y - Axis', fontsize=15)
    fig.colorbar(strasse, shrink=0.5, label='Height above the ground', aspect=15);


