.. jupyter-execute::

   #====================
   # 2d_pendulum_white_noise
   #====================
.. note::

   You can download this example as a Python script:
   :jupyter-download:script:`2d_pendulum_white_noise` or Jupyter notebook:
   :jupyter-download:notebook:`2d_pendulum_white_noise`.


.. jupyter-execute::


Basic 2D n particle pendulum, hanging from the ceiling. The point on the
ceiling may be fixed or move along the X and Y directions in a
predescribed manner. On each particle, white noise external forces are
applied in X and in Y direction.

sdeint is a library, which performs ITO (and Stratonovich) numerical
integration of stochastic differential equations. In the description it
says, it is still a ‘beta’ version. Playing around with it, it seems to
give ‘reasonable’ results.

.. jupyter-execute::

    import sympy.physics.mechanics as me
    import sympy as sm
    
    from scipy.optimize import fsolve
    import numpy as np
    
    import matplotlib.pyplot as plt
    %matplotlib inline
    from matplotlib import animation
    from IPython.display import HTML
    import matplotlib
    matplotlib.rcParams['animation.embed_limit'] = 2**128
    from matplotlib.patches import Rectangle
    
    import time
    
    import sdeint as sd

n is the number of particles. For simplicity, they all have the same
mass m, and the massless rods connecting them all have the same length
l.

Antrieb = True means the point on the ceiling moves in some predescribed
manner. I let it go around in a circle, other drives may be done.
Antrieb = False means the point on the ceiling is fixed.

All else is just standard sympy.physics.mechanics. - FX, FY hold the
external forces on each particle in X / Y direction. - reibung is the
friction in the joints. - kraft is the factor with which the external
forces may be multiplied to make them ‘stronger’. Same factor for all
particles. - rhs: to get the reaction forces one needs $ rhs = MM^{-1}
\* force $. This will be calculated numerically, hence these are just
‘substitutes’.

.. jupyter-execute::

    #----------------------------------------------------------------------------------------------------------
    n = 5                          #Number of particles.
    Antrieb = False                # if True the ceiling point moves in in pre described manner, else fixed
    #----------------------------------------------------------------------------------------------------------------------------------
    start = time.time()
        
    # Generalized coordinates and velocities
    q = []
    qd = []
    u = []
    rhs = []    # needed for the reaction force
    FX = []
    FY = []
    
    for i in range(1, n+1):                              #NOTE: data for particle i are at list place i-1
        q.append(me.dynamicsymbols('q'+str(i)))          # rotation of particle i relative to N
        u.append(me.dynamicsymbols('u'+str(i)))         
        qd.append(me.dynamicsymbols('q'+str(i), 1))
        rhs.append(me.dynamicsymbols('rhs'+str(i)))      # needed later for the reaction force
        FX.append(me.dynamicsymbols('F' + str(i)+ 'x'))  # external forces
        FY.append(me.dynamicsymbols('F' + str(i) + 'y')) #     dto.
        
    
    ux01, uy01, fx, fy = me.dynamicsymbols('ux01, uy01, fx, fy')          # reaction forces on suspension point
    aux = [ux01, uy01]
    
    g, l, m, reibung, kraft, t = sm.symbols('g, l, m, reibung, kraft, t') #kraft = intensity of white noise
    
    # of  course  any other 'reasonable' functions will also work
    if Antrieb == True:
        antriebx = l * sm.sin(sm.pi/2. * t)
        antriebxdt = antriebx.diff(t)
        antrieby = l * sm.cos(sm.pi/2. * t)
        antriebydt = antrieby.diff(t)
    else:
        antriebx = sm.S(0.)
        antriebxdt = sm.S(0)
        antrieby = sm.S(0.)
        antriebydt = sm.S(0.)
        
    #Reference frame, reference point
    N = me.ReferenceFrame('N')
    P0 = me.Point('P0')
    P0.set_vel(N, 0)
    
    P01 = P0.locatenew('P01', antriebx * N.x + antrieby * N.y)             # Base point on ceiling, 
    P01.set_vel(N, antriebxdt*N.x + antriebydt*N.y + ux01*N.x + uy01*N.y)  # ux01, uy01 for reaction forces
    
    # Frame points, particles for each Pendulum, numbered P01a, 1, 2, ..., n
    A = [N]
    P = [P01]
    particles = [me.Particle('P01a', P01, m)]                      # Particle fixed to ceiling
    
    for i in range(1, n+1):
        Ai = N.orientnew('A' + str(i), 'Axis', [q[i-1], N.z]) 
        Ai.set_ang_vel(N, u[i-1] * N.z)                      
        A.append(Ai)
        
        Pi = P[i-1].locatenew('P' + str(i), (l + (i-1)/n) * A[i].y)
        Pi.v2pt_theory(P[i-1], N, A[i])
        P.append(Pi)
        
        # Create a new particle of mass m at this point
        Pai = me.Particle('Pa' + str(i), Pi, m)
        particles.append(Pai)
    
    # kinematic equations
    kd = []
    for i in range(n):
        kd.append((qd[i] - u[i]))
    
    # forces
    FL = [(P01, - m*g*N.y + fx*N.x + fy*N.y)]  # point on ceiling different from the rest.
    for i in range(1, n+1):
        kraft_auf_punkt =(P[i], -m*g*N.y + kraft*(FX[i-1]*A[i].x + FY[i-1]*A[i].y))
        torque = (A[i], -reibung*u[i-1]*A[i].z)    
        FL.append(kraft_auf_punkt)
        FL.append(torque)
    
    # Kane's equations
    KM = me.KanesMethod(N, q_ind=q, u_ind=u, u_auxiliary=aux, kd_eqs=kd)
    (fr, frstar) = KM.kanes_equations(particles, FL)
    
    MM = KM.mass_matrix_full
    force = KM.forcing_full
    


A stochastic differential equation of the Ito type can have the random
forces only in a linearized manner, see any textbook on stochastic
differential equations. Hence the force vector must be linearized around
:math:`(FX, FY)`.

I think, this makes physical sense: if the force :math:`k` is small,
then :math:`f(x, k) =_{approx} f(x, 0) + d/dk(f(x, 0) * k`. Now here we
have: :math:`d/dt(x) = f(x, 0) + d/dkf(x, k) * k` or, as usually written
in stochastic textbooks $ dx_t = f(x, 0) \* dt + d/dkf(x, 0) \* dk_t$,
with :math:`dk_t` = white noise

:math:`forceito = force(FX=0, FY=0) + d/d(FX+FY)(force(FX, FY)`. Note,
that FX, FY are lists, and the diffentiation is done using the Jacobian.

In the code below, :math:`d/d(FX+FY)(force(FX, FY)` is called :math:`B`,
:math:`force(FX=0, FY=0)` is called :math:`force`, and :math:`forceito`
is not needed explicitly anywhere. Also some notational inconsitency on
my part: Kane’s formalism gives :math:`MM`, and :math:`force`, so that:
:math:`MM * d/dt(x) = force`. But in the integration we need $ d/dt(x) =
MM^{-1} \* (force + B) $

As :math:`FX, FY` are needed explicitly to get the rection forces at the
point on the ceiling, the mass matrix MM, and the force vectors have to
be enlarged properly.

.. jupyter-execute::

    # expand MM properly
    Hilfs = sm.Matrix.zeros(len(force), len(FX + FY))
    MM = sm.Matrix.hstack(MM, Hilfs)
    Hilfs = sm.Matrix.hstack(sm.Matrix.zeros(len(FX + FY), len(force)), sm.eye(len(FX + FY)))
    MM = sm.Matrix.vstack(MM, Hilfs)
    print('MM shape:', MM.shape)
    print('MM DS:', me.find_dynamicsymbols(MM))
    print('MM free symbols', MM.free_symbols)
    print('MM has {} operations'.format(np.sum(np.array([MM[i, j].count_ops(visual=False) for i in range(MM.shape[0]) for j in range(MM.shape[1])]))), '\n')
    
    # Linearize force around F = 0, so ITO's formalism may be applied
    B = force.jacobian(sm.Matrix([FX+FY]))
    #Extension as Fx, Fy are explizitly needed for the reaction forces on P01: fx, fy
    B = sm.Matrix.vstack(B, sm.Matrix(sm.eye(len(FX + FY))))
    print('B shape:', B.shape)
    print('B DS:', me.find_dynamicsymbols(B))
    print('B free symbols', B.free_symbols)
    print('B has {} operations'.format(np.sum(np.array([B[i, j].count_ops(visual=False) for i in range(B.shape[0]) for j in range(B.shape[1])]))), '\n')
    
    subs_dict = {i: sm.S(0.) for i in FX + FY}
    force = sm.Matrix.vstack(force, sm.Matrix([sm.S(0.) for i in range(len(FX+FY))])).subs(subs_dict) 
    print('force shape:', force.shape)
    print('force DS', me.find_dynamicsymbols(force))
    print('force has {} operations'.format(np.sum(np.array([force[i].count_ops(visual=False) for i in range(len(force))]))), '\n')


Set up the reaction forces. $ rhs = MM^{-1} \* force $ is needed for the
reaction forces. This will be calculated numerically later.

orte_x, orte_y are needed for the animation only.

I always find it interesting to look at the energy of the system.
Strange behaviour may point to mistakes in setting up Kane’s equations
of motion.

Lastly the sympy functions are converted into numpy functions, using
sm.lambdify(…)

.. jupyter-execute::

    # Reaction force
    subs_dict = {i.diff(t): rhs[j] for j, i in enumerate(u)}
    eingepraegt = KM.auxiliary_eqs.subs(subs_dict)                        
    print('eingepraegt shape:', eingepraegt.shape)
    print('eingepraegt DS:', me.find_dynamicsymbols(eingepraegt))
    print('eingepraegt free symbols', eingepraegt.free_symbols)
    print('eingepraegt has {} operations'.format(np.sum(np.array([eingepraegt[i].count_ops(visual=False) for i in range(eingepraegt.shape[0])]))), '\n')
    
    
    # needed for animation only
    orte_x = []
    orte_y = []
    
    for i in range(n+1):
        orte_x.append(me.dot(P[i].pos_from(P0), N.x))
        orte_y.append(me.dot(P[i].pos_from(P0), N.y))
    
    
    pot_energie = sum([m*g*me.dot(P[i].pos_from(P0), N.y) for i in range(n+1)])
    kin_energie = sum([particles[i].kinetic_energy(N).subs({i: 0. for i in aux}) for i in range(n+1)])
    
    #Lambdification
    qL = q + u + FX + FY
    pL = [g, m, l, reibung, kraft]
    F = [fx, fy]
    
    force_lam = sm.lambdify(qL + [t] + pL, force, cse=True)
    MM_lam = sm.lambdify(qL + [t] + pL, MM, cse = True)
    B_lam = sm.lambdify(qL + [t] + pL, B, cse=True)
    
    eingepraegt_lam = sm.lambdify(F + qL + [t] + pL + rhs, eingepraegt, cse=True)
    
    orte_x_lam = sm.lambdify(qL + [t] + pL, orte_x, cse=True)
    orte_y_lam = sm.lambdify(qL + [t] + pL, orte_y, cse=True)
    
    pot_lam = sm.lambdify(qL + [t] + pL, pot_energie, cse=True)
    kin_lam = sm.lambdify(qL + [t] + pL, kin_energie, cse=True)
    
    print("It took {:.3f} sec to establish Kane's formalism".format(time.time() - start))

For integrating this Ito type differential equation, I use the sdeint
library. It is the only one I know. It seems to be basically a ‘one man
enterprise’, not updated frequently, but at the moment it works. As to
the meaning of the parameters of itoint, consult the documentation of
sdeint.

sdeint is not optimized for speed, as the documentation says. Based on
my limited experience with it, I think, about 1,000 steps per second of
integration is about right. ( I asked the developer this question, but
never received a reply)

Input variables: - :math:`m1`: mass of each particle - :math:`l1`:
length of the massless rod - :math:`reibung1`: is the friction in the
joints - :math:`kraft1`: is the factor by which the ‘white noise’ is
multiplied

.. jupyter-execute::

    # Numerical integration
    start = time.time()
    
    # Input variables
    #==================================================================================================
    
    intervall = 3.                       # integration runs from 0. to intervall
    # number of evaluations of sdeint. One should have at least 500 evaluations per time unit.
    schritte = 3000
    m1 = 1.                             # mass of each particle
    l1 = 1.                             # distance from one particle to the next
    reibung1 = 0.                       # friction 
    kraft1 = 3.                         # 'strength' of forces
    #==================================================================================================
    
    times = np.linspace(0, intervall, schritte)
    
    pL_vals = [9.8, m1, l1, reibung1, kraft1]                
    
    # for simplicity, the initial conditions are always the same:
    # pendulum is hanging straight down, external forces are 0. at t = 0.
    y0 = [np.pi for i in range(n)] + [0. for i in range(n, 2*n)] + [0.for i in range(len(FX + FY))]                      #Anfangsbed.
    
    # the 'meaning' of f, G may be learned from sdeint site.
    def f(x, t):
        A = np.linalg.solve(MM_lam(*x, t, *pL_vals), force_lam(*x, t, *pL_vals)).reshape(len(force))
        return A 
    
    def G(x, t):
        B = np.linalg.solve(MM_lam(*x, t, *pL_vals), B_lam(*x, t, *pL_vals)) 
        return B
    
    resultat = sd.itoint(f, G, y0, times)
    
    print('it took {:.3f} sec to integrate {} steps'.format(time.time() - start, schritte))
    print('resultat.shape: ', resultat.shape)
    
        
    
    #====================================================================================


Print the angular speeds. I guess, all it shows is that they look
‘random’, as they should.

.. jupyter-execute::

    #==================================================================================
    # print generalized angular velocities
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, j in zip(range(n, 2*n), [i for i in range(1, n+1)]):
        ax.plot(times, resultat[:, i], label='Angular velocity of node {}'.format(j))
        ax.set_title('Generalized angular velocities')
    ax.legend();


Print the externally applied forces. They should look like Brownian
motion, integrated white noise, and they seem to look o.k.

.. jupyter-execute::

    start = time.time()
    
    #Print externally applied forces 
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i, j in zip(range(2*n, 4*n), FX + FY):
        ax.plot(times, kraft1*resultat[:, i], label=' {}'.format(j))
        ax.set_title('Integrated externally applied forces, that is Brownian Motion')
    ax.legend();


Calculate the reaction forces on the suspension point on the ceiling. -
numerically calculate rhs - numerically solve eingepraegt for fx, fy -
plot them

.. jupyter-execute::

    # print reaction forces on P01
    # 1. numerically calculate the RHS, as it is needed for the reaction forces.
    RHS = np.zeros((schritte, len(force)))
    for i in range(schritte):
        zeit = times[i]
        RHS[i, :] = np.linalg.solve(MM_lam(*[resultat[i, j]for j in range(resultat.shape[1])], zeit, *pL_vals), 
                                   force_lam(*[resultat[i, j] for j in range(resultat.shape[1])], zeit, *pL_vals)).reshape(resultat.shape[1])
    
    print('RHS shape', RHS.shape)
    
    # 2. calculate reaction forces numerically: solve eingepraegt for fx, fy
    kraft_x = np.zeros(schritte)
    kraft_y = np.zeros(schritte)
    
    def func (x, *args):
    # just needed to 'convert the arguments' properly
        return eingepraegt_lam(*x, *args).reshape(2)
    
    x0 = tuple([1., 1.])   #initial guess
    
    for i in range(schritte):
        zeit = times[i]
        y0 = [resultat[i, j] for j in range(resultat.shape[1])]
        rhs = [RHS[i, j] for j in range(n, 2*n)]
        pL_vals = list(pL_vals)
        
        args = tuple(y0 + [zeit] + pL_vals + rhs)
        A = fsolve(func, x0, args=args).reshape(2)
        x0 = tuple(A)      # new guess, should speed up convergence.
        kraft_x[i] = A[0]
        kraft_y[i] = A[1]
    
    
    
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, kraft_x, label='Reaction force on P0 in X direction')
    ax.plot(times, kraft_y, label='Reaction force on P0 in Y direction')
    ax.set_title('Reaction forces on P01')
    ax.legend();


Plot the energies of the system. As is to be expected, they increase /
and decrease randomly. As is to be expected, the potential energy is
smoother than the kinetic energy.

.. jupyter-execute::

    pot_np = np.empty(schritte)
    kin_np = np.empty(schritte)
    total_np = np.empty(schritte)
    
    for i in range(schritte):
        zeit = times[i]
        pot_np[i] = pot_lam(*[resultat[i, j] for j in range(resultat.shape[1])], zeit, *pL_vals)
        kin_np[i] = kin_lam(*[resultat[i, j] for j in range(resultat.shape[1])], zeit, *pL_vals)
        total_np[i] = pot_np[i] + kin_np[i]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, pot_np, label='Potential Energy')
    ax.plot(times, kin_np, label='Kinetic Energy')
    ax.plot(times, total_np, label='Total Energy')
    
    ax.set_title('Energy of the System of {} bodies'.format(n+1))
    ax.legend();
    
    
    print('it took {:.3f} sec to plot the graphs'.format(time.time() - start))

Create an animation.

The last comand HTML(..) is needed on my iPad to show an animation. I do
not know, whether needed on other machines.

.. jupyter-execute::

    
    # This is only to reduce the number of framesfor the animation, else it takes forever on my iPad
    # with HTML
    if schritte > 1000:
        N = int(schritte / 200)  #to make the animation faster
    else:
        N = 1
        
    
    x_coords = np.zeros((int(schritte/N)+1, n+1))
    y_coords = np.zeros((int(schritte/N)+1, n+1))
    
    k = -1
    for j in range(len(times)):
        zeit= times[j]
        if j % N == 0:
            k += 1
            x_coords[k] = orte_x_lam(*[resultat[j, i] for i in range(resultat.shape[1])], zeit, *pL_vals )
            y_coords[k] = orte_y_lam(*[resultat[j, i] for i in range(resultat.shape[1])], zeit, *pL_vals )
        
    def animate_pendulum(times, x, y):
        
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.axis('on')
        lim1 = max([np.abs(y_coords[i, n]) for i in range(int(schritte/N))]) + l1/5.
        lim2 = max([np.abs(y_coords[i, n]) for i in range(int(schritte/N))]) + l1/5.
        lim = max(lim1, lim2)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim))
    
        line, = ax.plot([], [], 'o-', lw=0.5)
        
        cart_width = 0.2 * lim
        cart_height = 0.1 * lim
        rect = Rectangle([(x[0, 0]) - cart_width/2, y[0, 0] - cart_height/2], 
                         cart_width, cart_height, fill=True, color='red', ec='black')
        ax.add_patch(rect)
        
        def init():
            line.set_data([], [])
            return line,
    
        def animate(i):
            rect.set_xy((x[i, 0] - cart_width/2., y[i, 0] - cart_height/2.))
            line.set_data(x[i], y[i])
            zeit = i * N * intervall / schritte 
            ax.set_title('stochastically exited pendulum at time {:.2f} sec'.format(zeit), fontsize=15)
            return line,
    
        anim = animation.FuncAnimation(fig, animate, frames=int(schritte/N),
                                       interval=1000*times.max() / (4*N) ,
                                       blit=True, init_func=init)
        plt.close(fig)
        return anim
    
    anim = animate_pendulum(times, x_coords, y_coords)
    HTML(anim.to_jshtml())
    


