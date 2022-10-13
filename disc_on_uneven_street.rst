
... note::

You can download this example as a Python script:
:jupyter-download:script: disc_on_uneven_stree or Jupyter notebook:
:jupyter-download:notebook: disc_on_uneven_street.

.. jupyter_execute::

    import sympy as sm
    import sympy.physics.mechanics as me
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve, minimize
    
    from matplotlib import animation
    import matplotlib
    from matplotlib import patches
    import matplotlib.pyplot as plt
    from IPython.display import HTML
    matplotlib.rcParams['animation.embed_limit'] = 2**128
    %matplotlib inline
    
    import time

.. code:: ipython3

    '''
    A homogeneous disc with radius r and mass m is running on an uneven 'street' without sliding. 
    The disc is not allowed to jump, hence I can calculate the reaction forces needed to hold it on the
    street at all times
    
    The 'overall' shape of the street is modelled as a parabola (called strassen_form), the unevenness is 
    modelled as a sum of sin functions (called strasse) with each term having a smaller amplitude and higher 
    frequency as the previous one.
    the 'street itself' is the sum of strasse and strassen_form.
    
    NOTE: The disc must always have only have one contact point, hence this inequality 
    must hold: r < abs( (1 + strasse.diff(t)**2)**3/2 / strasse.diff(t, 2) ), 
    formula for bevel radius (Krümmungsradius) from Wikipedia.
    
    Using simplify() on intermediate results (which are small enough so it ends in good time!) 
    reduce the number of operations substantially.
    
    The total energy is not constant. I think, these are numerical inaccuracies, as this becomes worse 
    if the radius of the disc is close to the 'bevel radius'.
    method = 'Radau' in solve_ivp seems to be best.
    '''
    
    start = time.time()
    
    # q1 is the angle of the disc, x is the horizontal position of contact point CP
    q1, x = me.dynamicsymbols('q1 x')  
    u1, ux = me.dynamicsymbols('u1 ux')
    
    
    auxx, auxy, fx, fy = me.dynamicsymbols('auxx auxy fx fy')   # for the reaction forces at contact point CP
    
    m, g, r, iZZ = sm.symbols('m, g, r, iZZ')
    amplitude, frequenz, reibung, t = sm.symbols('amplitude frequenz reibung t') # reibung = friction in German
    
    
    N = me.ReferenceFrame('N')         # fixed inertial frame
    A1 = me.ReferenceFrame('A1')       # street. A1.x will be perpendicular to the derivative of the street at CP
    A2 = me.ReferenceFrame('A2')       # fixed to the disc
    
    P0 = me.Point('P0')                # fixed point
    CP = me.Point('CP')                # contact point of wheel and street
    Dmc = me.Point('Dmc')              # geometric center of the disc
    
    # make the uneven street.
    #============================================
    rumpel = 3    # the higher the number the more 'uneven the street'
    #============================================
    strasse = sum([amplitude/j * sm.sin(j*frequenz * x) for j in range(1, rumpel)])
    strassen_form = (frequenz/1.5 * x)**2
    gesamt = strassen_form + strasse
    #===========================================
    
    '''
    Relationship of x(t) to q1:
    x(t) is the horizontal position of the contact point.
    The arc length of a function f(k) from 0 to x is integral[0, x](sqrt(1 + f(k).diff(t)**2)dk) [I found this
    in the internet]
    If the disc turns by an angle q1, the contact point moves a distance q1 * r , hence
    r * (-q(t)) = arc length of f(k)[0, x(t)]  = integral[0, x](sqrt(1 + f(k).diff(k)**2)dk) [the minus sign
    is 'dicated' by the right hand rule].  Hence:
    d/dt(r * q) = r * d/dt(q(t)) =  d/dt(arc length(..)) = sqrt(1 + strasse.diff(x)**2) * d/dt(x). 
    hence:
    d/dt(x) = -r * u / sqrt(1 + strasse.diff(x)**2) is the sought after differential equation to determine x(q(t)).
    '''
    rhs3 = (-u1 * r / sm.sqrt(1. + gesamt.diff(x)**2)).simplify()
    
    
    '''
    this is to determine the maximum radius of the wheel so it have only one contact point.
    the radius of the wheel must be smaller than the smallest bevel radius (Schmiegekreis) of the
    function of the road. 
    '''
    r_max = (sm.S(1.) + (gesamt.diff(x))**2 )**sm.S(3/2)/gesamt.diff(x, 2)
    
    '''
    The vector perpendicular to the strasse is -(gesamt.diff(x), - 1). The leading minus sign, because directed
    'inward'. It points from the contact point CP to the geometric center of the discDmc
    '''
    vector = (-(gesamt.diff(x)*N.x - N.y)).simplify()
    A2.orient_axis(N, q1, N.z)
    A2.set_ang_vel(N, u1 * N.z)
    
    CP.set_pos(P0, x*N.x + gesamt*N.y)    # location of contact point
    # CP has no real velocity, as it is 'part of the street' 
    # auxx, auxy are virtual speeds to determine the reaction forces on CP
    CP.set_vel(N, auxx*N.x + auxy*N.y)
    CP_pos = [me.dot(CP.pos_from(P0), uv) for uv in (N.x, N.y)]
    
    
    #The center of the wheel is at distance r from CP, perpendicular to the surface of the street.
    dir_Dmc = (vector.normalize()).simplify()
    Dmc.set_pos(CP, r * dir_Dmc)
    
    # for an 'instant' CP is fixed in A2, as CP has zero speed, being part of the street, too.
    # Hence I can use v2pt_theory as below to get the speed of Dmc.
    Dmc.v2pt_theory(CP, N, A2)
    Dmc_pos = [me.dot(Dmc.pos_from(P0), uv) for uv in (N.x, N.y)]
    
    I = me.inertia(A2, 0., 0., iZZ)                                              
    Body = me.RigidBody('Body', Dmc, A2, m, (I, Dmc))                               
    BODY = [Body]
    
    '''
    A necessary, but by no means sufficient condition for the correctness of the equations of motion is that,
    absent any friction, the total energy be constant.
    Hence I like to look at this.
    '''
    kin_energie = Body.kinetic_energy(N).subs({auxx: 0., auxy:  0.})
    pos_energie = m * g * me.dot(Dmc.pos_from(P0), N.y)
    
    # Setting up Kane's formalism.
    #=======================================================================
    FL = [(Dmc, -m*g*N.y), (CP, fx*N.x + fy*N.y), (A2, -reibung*u1*A2.z)]
    kd = [u1 - q1.diff(t)]  # kinematic equations
    
    q = [q1]
    u = [u1]
    aux = [auxx, auxy]
    
    KM = me.KanesMethod(N, q_ind=q, u_ind=u, kd_eqs=kd, u_auxiliary=aux)
    (fr, frstar) = KM.kanes_equations(BODY, FL)
    MM = KM.mass_matrix_full
    force = KM.forcing_full
    '''
    rhs is needed for the reaction forces. Here it is small enough to use it symbolically. If it becomes
    large, it is better to calculate it numerically.
    '''
    rhs = KM.rhs().subs({sm.Derivative(x, t): rhs3})
    print('rhs DS', me.find_dynamicsymbols(rhs))
    print('rhs free symbols', rhs.free_symbols)
    print('rhs has {} operations'.format(sum([rhs[i].count_ops(visual=False) for i in range(len(rhs))])), '\n')
    
    
    # Reaction forces
    eingepraegt = KM.auxiliary_eqs.subs({sm.Derivative(u1, t): rhs[1], sm.Derivative(x, t): rhs3})
    print('eingepraegt DS', me.find_dynamicsymbols(eingepraegt))
    print('eingepraegt free symbols', eingepraegt.free_symbols)
    print('eingepraegt has {} operations'.format(sum([eingepraegt[i].count_ops(visual=False) for i in range(len(eingepraegt))])), '\n')
    
    # Add rhs3 at the bottom of force, to get d/dt(x) = rhs3. This is to numerically integrate x(t)
    force = sm.Matrix.vstack(force, sm.Matrix([rhs3])).subs({sm.Derivative(x, t): rhs3})
    print('force DS', me.find_dynamicsymbols(force))
    print('force free symbols', force.free_symbols)
    print('force has {} operations'.format(sum([force[i].count_ops(visual=False) for i in range(len(force))])), '\n')
    
    # Enlarge MM properly
    MM = sm.Matrix.hstack(MM, sm.Matrix([0., 0.]))
    MM = sm.Matrix.vstack(MM, sm.Matrix([0., 0., 1.]).T)
    print('MM DS', me.find_dynamicsymbols(MM))
    print('MM free symbols', MM.free_symbols)
    print('MM has {} operations'.format(sum([MM[i, j].count_ops(visual=False) for i in range(MM.shape[0]) for j in range(MM.shape[1])])), '\n')
    
    
    
    # Lambdification. Turning symbolic expressions into numpy functions.
    pL = [m, g, r, iZZ, amplitude, frequenz, reibung]
    qL = q + u + [x]
    F = [fx, fy]
    
    MM_lam = sm.lambdify(qL + pL, MM, cse=True)
    force_lam = sm.lambdify(qL + pL, force, cse=True)
    
    CP_pos_lam = sm.lambdify(qL + pL, CP_pos, cse=True)
    Dmc_pos_lam = sm.lambdify(qL + pL, Dmc_pos, cse=True)
    
    # will be solved for F numerically later. Much too large to be solved symbollically.
    eingepraegt_lam = sm.lambdify(F + qL + pL, eingepraegt, cse=True) 
    
    #this is needed to plot the shape of the street
    strasse_lam = sm.lambdify([x] + pL,  gesamt, cse=True)
    
    kin_lam = sm.lambdify(qL + pL, kin_energie, cse=True)
    pos_lam = sm.lambdify(qL + pL, pos_energie, cse=True)
    
    r_max_lam = sm.lambdify([x] + pL, r_max,cse=True)
    
    print('it took {:.3f} sec to establish Kanes equations'.format(time.time() - start))


.. code:: ipython3

    # Integrate numerically
    
    start = time.time()
    
    # Input parameters 
    #==========================================================
    mm = 1.
    rr = 4.
    amplitude = 1.
    frequenz = 0.25     # the smaller this number, the more 'even' the street   
    reibung = 0.        # Friction
    intervall = 25.     # time inverval of integration is [0., intervall]
    schritte = 500      # Where the results of solve_ivp will be given, see description of solve_ivp.
    
    q0 = 0.             # starting angle. As the disc is symmetric about this angle, it plays no real role
    u0 = 8.             # starting angular velocity of disc.
    x0 = 0.             # Starting X position of disc. 
    #==========================================================
    
    
    iZZe = 1/2 * mm * rr**2
    pL_vals = [mm, 9.8, rr, iZZe, amplitude, frequenz, reibung]
    y0 = [q0, u0, x0]
    print('Arguments')
    print('[m, g, r, iXX, iYY, iZZ, iXY, iXZ, iYZ, amplitude, frequenz, reibung]')
    print(pL_vals, '\n')
    print('[q0, u0, x0]')
    print(y0, '\n')
    
    startwert = y0[2]   # just needed for the plots below
    startomega = y0[1]  #  dto.
    
    #find the largest admissible r, given strasse, amplitude, frequenz
    def func(x, args):
    # just needed to get the arguments matching for minimize
        return np.abs(r_max_lam(x, *args))
    
    x0 = 0.1            # initial guess
    minimal = minimize(func, x0, pL_vals)
    
    if pL_vals[2] < (x := minimal.get('fun')):
        print('selected radius = {} is less than maximally admissible radius = {:.2f}, hence o.k.'.format(pL_vals[2], x), '\n')
    else:
        print('selected radius {} is larger than admissible radius {:.2f}, hence NOT o.k.'.format(pL_vals[2], x), '\n')
        
        
    times = np.linspace(0, intervall, schritte)
                            
    def gradient(t, y, args):
        vals = np.concatenate((y, args))
        sol = np.linalg.solve(MM_lam(*vals), force_lam(*vals))
        return np.array(sol).T[0]
    
    # method = 'Radau' seems to work best here.
    resultat1 = solve_ivp(gradient, (0., float(intervall)), y0, t_eval=times, args=(pL_vals,), 
                method='Radau', atol=1.e-12, rtol=1.e-9)
    resultat = resultat1.y.T
    event_dict = {-1: 'Integration failed. Do not run the plot', 0: 'Integration finished successfully', 1: 'some termination event'}
    print(event_dict[resultat1.status])
    print('resultat shape', resultat.shape, '\n')
    
    print("To numerically integrate an intervall of {} sec the routine cycled {} times and it took {:.5f} sec ".format(intervall, resultat1.nfev, time.time() - start))


.. code:: ipython3

    # plot results
    
    start = time.time()
    
    Dmc_X = np.empty(schritte)
    Dmc_Y =np.empty(schritte)
    for i in range(schritte):
        Dmc_X[i], Dmc_Y[i] = Dmc_pos_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals)
    
    
    fig, ax = plt.subplots(figsize=(15, 5))
    for i, j in zip(range((resultat.shape[1])), ('rotational angle', 'rotational speed', 'displacement')):
        ax.plot(times, resultat[:, i], label=j)
    ax.set_title('Coordinates')
    ax.legend();
    
    #calculate implied forces numerically
    def func (x, *args):
    # just serves to make the arguments compatible between fsolve and eingepraegt_lam
        return eingepraegt_lam(*x, *args).reshape(len(F))
    
    kraft = np.zeros((schritte, len(F)))
    x0 = tuple([1. for i in range(len(F))])   # initial guess
    for i in range(schritte):
        y00 = [resultat[i, j] for j in range(resultat.shape[1])]
        args = tuple((y00 + pL_vals))
        A = fsolve(func, x0, args=args).reshape(len(F)) # numerically find fx, fy
        x0 = tuple(A)      # updated initial guess, should improve convergence
        kraft[i] = A        
            
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(times, kraft[:, 0], label = 'Fx')
    ax.plot(times, kraft[:, 1], label = 'Fy')
    ax.set_title('Reaction forces on contact point')
    ax.legend();
    
    
    fig, ax = plt.subplots(figsize=(15, 5))
    
    # plot the street, and the extremes, of the position of the disc.
    links = np.min(resultat[:, 2])
    rechts = np.max(resultat[:, 2])
    ruhe = np.mean([resultat[-30::, 2]])    # get approx. rest position of wheel
    maximal = max(np.abs(links), np.max(rechts))
    times1 = np.linspace(-maximal-5, maximal+5, schritte)
    ax.plot(times1, strasse_lam(times1, *pL_vals)  , label='Strasse')
    if pL_vals[-1] != 0.:
        ax.axvline(ruhe,ls = '--', color='red', label='approx. fimal pos. of wheel')
    ax.axvline(links,ls = '--', color='green', label='leftmost pos. of wheel')
    ax.axvline(rechts,ls = '--', color='black', label='rightmost pos. of wheel');
    ax.axvline(startwert, ls='--', color='orange', label='starting position of wheel')
    if startomega > 0.:
        richtung = 'left'
    else:
        richtung = 'right'
    text = 'Wheel has speed ' + str(np.abs(startomega)) + ' units to the ' + richtung
    plt.title(text)
    ax.legend();
    #===========
    
    kin_np = np.empty(schritte)
    pos_np = np.empty(schritte)
    total_np = np.empty(schritte)
    
    for i in range(schritte):
        kin_np[i] = kin_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals)
        pos_np[i] = pos_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals)
        total_np[i] = kin_np[i] + pos_np[i]
    
    if pL_vals[-1] == 0.:
        print('Max deviation from constant of total energy is {:.4f} %'.format((max(total_np) - min(total_np))/max(total_np) * 100.))
    
    fig, ax = plt. subplots(figsize=(15, 5))
    ax.plot(times, kin_np, label='kinetic energy')
    ax.plot(times, pos_np, label='pos energy')
    ax.plot(times, total_np, label='total energy')
    ax.set_title('Energy of the disc')
    ax. legend();
    print('it took {:.3f} sec to calculate the forces and plot the graphs'.format(time.time() - start))


.. code:: ipython3

    #Animation
    # Location of the center of the disc
    Dmcx = np.empty(schritte)
    Dmcy =np.empty(schritte)
    for i in range(schritte):
        Dmcx[i], Dmcy[i] = Dmc_pos_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals)
    
    # needed to give the picture the right size.
    xmin = min([resultat[i, 2] for i in range(schritte)])
    xmax = max([resultat[i, 2] for i in range(schritte)])
    
    ymin = min([strasse_lam(resultat[i, 2], *pL_vals) for i in range(schritte)]) 
    ymax = max([strasse_lam(resultat[i, 2], *pL_vals) for i in range(schritte)]) 
    
    # Data to draw the uneven street
    cc = rr
    strassex = np.linspace(xmin - 3*cc, xmax + 3.*cc, schritte)
    strassey = [strasse_lam(strassex[i], *pL_vals) for i in range(len(strassex))]
    
    # Data for the dashed lines, which mark the contact point of the ellipse with the street
    geradex = np.linspace(xmin - 3.*cc, xmax + 3.*cc, schritte)
    geradey = np.linspace(ymin - 3.*cc, ymax + 3.*cc, schritte)
    
    vertikal = np.empty(schritte)
    horizontal =np.empty(schritte)
    vertikal_store = []
    horizontal_store = []
    pointx = []
    pointy = []
    
    for i in range(schritte):
        vertikal = np.array([resultat[i, 2] for k in range(schritte)])
        vertikal_store.append(vertikal)
        horizontal = np.array([strasse_lam(resultat[i, 2], *pL_vals) for k in range(schritte)])
        horizontal_store.append(horizontal)
        pointx.append(Dmcx[i] + rr*np.cos(resultat[i, 0]))  #data for the point of the disc.
        pointy.append(Dmcy[i] + rr*np.sin(resultat[i, 0]))  #      dto.
        
        
        
    def animate_pendulum(times, x1, y1):
        
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'aspect': 'equal'})
        
        ax.axis('on')
        ax.set_xlim(xmin - 3.*cc, xmax + 3.*cc)
        ax.set_ylim(ymin - 3.*cc, ymax + 3.*cc)
        ax.plot(strassex, strassey)
    
    
        line1, = ax.plot([], [], 'o-', lw=0.5)
        line2, = ax.plot([], [], linestyle = '--')
        line3, = ax.plot([], [], linestyle = '--')
        line4, = ax.plot([], [], 'bo', markersize=10) # the dot on the disc, to show it is rotating
        
        elli = patches.Circle((x1[0], y1[0]), radius = rr, fill=True, color='red', ec='black')
        ax.add_patch(elli)
    
        def animate(i):
            
            ax.set_title('running time {:.2f} sec'.format(times[i]), fontsize=15)
            
            elli.set_center((x1[i], y1[i]))
            elli.set_height(2.*rr)
            elli.set_width(2.*rr)
            elli.set_angle(np.rad2deg(resultat[i, 0]))
                           
            line1.set_data(x1[i], y1[i])                  # center of the ellipse
            line2.set_data(geradex, horizontal_store[i])  # dashed line to mark the contact point
            line3.set_data(vertikal_store[i], geradey)    #            dto. 
            line4.set_data(pointx[i], pointy[i])
            return line1, line2, line3, line4,
    
        anim = animation.FuncAnimation(fig, animate, frames=len(times),
                                       interval=1000*max(times) / len(times),
                                       blit=True)
        plt.close(fig)
        return anim
    
    anim = animate_pendulum(times, Dmcx, Dmcy)
    #HTML(anim.to_jshtml())    # needed, when run on an iPad, I know no other way to do it. It is SLOW!
    


