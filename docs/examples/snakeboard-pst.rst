.. jupyter-execute::

#===============
# snakeboard-pst
#===============

.. note::

    You can download this example as a Python script:
    :jupyter-download:script:'snakeboard-pst' or Jupyter notebook:
    :jupyter-download:notebook:'snakeboard-pst'.

.. jupyter-execute::
    
    import sympy as sm
    import sympy.physics.mechanics as me
    import time
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve
    import matplotlib.pyplot as plt
    import matplotlib as mp
    from mpl_toolkits.mplot3d import axes3d
    %matplotlib inline
    from IPython.display import HTML
    mp.rcParams['animation.embed_limit'] = 2**128 # I need this for large animations on my iPad
    from matplotlib import animation

This is basically the snakeboard of Jason Moore’s lecture:
https://moorepants.github.io/learn-multibody-dynamics/motion.html
https://moorepants.github.io/learn-multibody-dynamics/_images/motion-snakeboard.svg

All I changed is this:

-  changed the variables to some extend, to use parts on an earlier
   program
-  added wheels, on which particles may be attached
-  changed the direction of the frames, so gravity points in the negtive
   N.y directon
-  used sympy.physics.mechanics to get the equations of motion with
   Kane’s method

This shows a sketch of my snake board:
https://postimg.cc/gallery/W3nJxw1 Click on it twice to remove the
advertizement.

===============================================================================================================

Define the sympy symbols and the dynamic symbols needed for the
calculation - :math:`q_1, q_2, q_3, q_4, q_5` are the coordinates
exactly like in JM’s example - :math:`u_1,... u_5` are the speeds -
:math:`q_{vL}, q_{vR}, q_{hL}, q_{hR}` are the angles of the wheels -
:math:`u_{vL}...` their speeds

-  the rest are constants needed to describe the physical properties of
   the system:

   -  reibung is the friction in the bearings of the wheels
   -  :math:`r, l_1, l_2, m, m_o, m_a, m_v` radius of the wheels, length
      of axles, length of ‘connecting’ axle, mass of wheel, mass of
      particles on the wheels, mass of axles, mass of connecting axle
   -  :math:`i_{XX}, i_{YY}, i_{ZZ}` are the moments of inertia of the
      wheel
   -  :math:`i_{XXa}, i_{YYa}, i_{ZZa}` are the moments of inertia of
      the axles
   -  :math:`i_{XXV}, i_{YYV}, i_{ZZV}` are the moments of inertia of
      the ‘connecting’ axle
   -  $ :raw-latex:`\alpha`\ *{vL}, :raw-latex:`\beta`*\ {vL},
      :raw-latex:`\alpha`\ *{vL}, :raw-latex:`\beta`*\ {vR}$, etc.
      determine the location of the observer w.r.t. the center of the
      wheel.

.. jupyter-execute::

    start = time.time()
    start1 = time.time()
    q1, q2, q3, q4, q5 = me.dynamicsymbols('q1, q2, q3, q4, q5') 
    u1, u2, u3, u4, u5 = me.dynamicsymbols('u1, u2, u3, u4, u5')
    qvL, qvR, qhL, qhR = me.dynamicsymbols('qvL, qvR, qhL, qhR')
    uvL, uvR, uhL, uhR = me.dynamicsymbols('uvL, uvR, uhL, uhR')
    
    mv, mh, mm = sm.symbols('mv, mh, mm', cls=me.Point)     # midpoints on the axles
            
    reibung, alphavL, betavL, alphavR, betavR = sm.symbols('reibung, alphavL, betavL, alphavR, betavR')
    alphahL, betahL, alphahR, betahR = sm.symbols('alphahL, betahL, alphahR, betahR')
    m, mo, g, r, l1, l2, iXX, iYY, iZZ, t = sm.symbols('m, mo, g, r, l1, l2, iXX, iYY, iZZ, t')
    ma, mV, iXXa, iYYa, iZZa, iXXV, iYYV, iZZV = sm.symbols('ma, mV, iXXa, iYYa, iZZa, iXXV, iYYV, iZZV')

Define the frames and the points as follows: - N = inertial frame -
:math:`A_{vX}` = frame attached to the front axle - :math:`A_V` = frame
attached to the link connecting the two axles - :math:`A_{hX}` = frame
attached to the rear axle - :math:`A_{vL}` = frame attached to the left
front wheel - :math:`A_{vR}` = frame attached to the right front wheel -
:math:`A_{hL}` = frame attached to the left back wheel - :math:`A_{hR}`
= frame attached to the right back wheel

-  :math:`P_0` = reference point, fixed in N

-  :math:`Dmc_{vL}` = center of the front left wheel

-  :math:`Dmc_{vR}` = center of the front right wheel

-  :math:`Dmc_{hL}` = center of the back left wheel

-  :math:`Dmc_{hR}` = center of the back right wheel

-  :math:`o_{vL}, o_{vR}` = particles anywhwere within the wheels

-  :math:`o_{hL}, o_{hR}` = particles anywhwere within the wheels

.. jupyter-execute::

    N, AvX, AvL, AvR, AV, AhX, AhL, AhR = sm.symbols('N, AvX, AvL, AvR, AV, AhX, AhL, AhR', cls = me.ReferenceFrame)         
    P0, DmcvL, DmcvR, DmchL, DmchR, ovL, ovR, ohL, ohR = sm.symbols('P0, DmcvL, DmcvR, DmchL, DmchR ovL, ovR, ohL, ohR', cls = me.Point)
    P0.set_vel(N, 0)
    
    AV.orient_axis(N, q3, N.y)   # frame of the axle connecting the two axles with wheels
    
    AvX.orient_axis(AV, q4, AV.y) # frame of front axle
    AhX.orient_axis(AV, q5, AV.y) # frame of rear axle
    
    AvL.orient_axis(AvX, qvL, AvX.x) # frames of the wheels
    AvR.orient_axis(AvX, qvR, AvX.x)
    AhL.orient_axis(AhX, qhL, AhX.x)
    AhR.orient_axis(AhX, qhR, AhX.x)

Define the various points, and their speeds. Here **v2pt_theory** is
probably the ideal method to get the speeds.

.. jupyter-execute::

    mm.set_pos(P0, q1*N.x + q2*N.z)     # center point of the axle bearing the two axles with wheels
    mm.set_vel(N, u1*N.x + u2*N.z)
    mv.set_pos(mm, -l2/2.*AV.z)         # center point of the front axle
    mv.v2pt_theory(mm, N, AV)
    mh.set_pos(mm, l2/2.*AV.z)          # center point of the rear axle
    mh.v2pt_theory(mm, N, AV)
    
    DmcvL.set_pos(mv, -l1/2.*AvX.x)      # center of mass of left front wheel
    DmcvL.v2pt_theory(mv, N, AvX)
    DmcvR.set_pos(mv, l1/2.*AvX.x)
    DmcvR.v2pt_theory(mv, N, AvX)
    DmchL.set_pos(mh, -l1/2.*AhX.x)
    DmchL.v2pt_theory(mh, N, AhX)
    DmchR.set_pos(mh, l1/2.*AhX.x)
    DmchR.v2pt_theory(mh, N, AhX)
    
    ovL.set_pos(DmcvL, r*alphavL * AvL.y + r*betavL * AvL.z)  # location of the particle on the front left wheel
    ovL.v2pt_theory(DmcvL, N, AvL)
    ovR.set_pos(DmcvR, r*alphavR * AvR.y + r*betavR * AvR.z)
    ovR.v2pt_theory(DmcvR, N, AvR)
    ohL.set_pos(DmchL, r*alphahL * AhL.y + r*betahL * AhL.z)
    ohL.v2pt_theory(DmchL, N, AhL)
    ohR.set_pos(DmchR, r*alphahR * AhR.y + r*betahR * AhR.z)
    ohR.v2pt_theory(DmchR, N, AhR);

Find the various **speed constraints** imposed by - :math:`m_v, m_h` are
not allowed to move in the :math:`A_{vX}`.x / :math:`A_{hX}`.x
directions respectively. This gives two speed constraints for the five
generalized speed :math:`u_1, ..u_5`. Like JM in his lecture, I solve
for :math:`u_1, u_2` making them the dependent generalized coordinates.
Here I use *sm.solve* to get the solution dictionary, for larger
systems, JM’s method using the *Jacobian* is much faster! - The
additional dependent speeds, :math:`u_{vL}, u_{vR}, u_{hR}` are needed
to describe the rotations of the wheels. Of course, they must rotate
with the speed determined by the speeds of the centers of the wheels. I
simply use the formula: *speed of point X on body B = angular velocity
of B*\ **x**\ *distance of X from a point (momentarily) at rest in B*,
where **x** indicates the cross product. To get the angular velocities
:math:`u_{vL}`, etc, I solve the equation above for the angular
velocity. I found this in the internet:
https://math.stackexchange.com/questions/2195047/solve-the-vector-cross-product-equation

.. jupyter-execute::

    # the speed of mv in AvX.x direction must be zero, dto for mh
    vel_mv_x = me.dot(mv.vel(N), AvX.x)
    vel_mh_x = me.dot(mh.vel(N), AhX.x)
    vel_m_dict = sm.solve((vel_mv_x, vel_mh_x), (u1, u2))
    
    # from the internet mentioned above
    uvL11 = (-1./r**2 * DmcvL.vel(N).cross(r*AvX.y)).subs(vel_m_dict)
    uvR11 = (-1./r**2 * DmcvR.vel(N).cross(r*AvX.y)).subs(vel_m_dict)
    uhL11 = (-1./r**2 * DmchL.vel(N).cross(r*AhX.y)).subs(vel_m_dict)
    uhR11 = (-1./r**2 * DmchR.vel(N).cross(r*AhX.y)).subs(vel_m_dict)
    '''
    #check, that these rotational speeds only have a component in AvX.x / AhX.x direction
    for i, j in zip((uvL11, uvR11), ('uvL11', 'uvR11')):
        test = i.express(AvX).subs(vel_m_dict).simplify()
        print(j +  ' only components in AvX.x direction? ', sm.Abs(me.dot(test, AvX.y)) + sm.Abs(me.dot(test, AvX.z)) == 0)
    for i, j in zip((uhL11, uhR11), ('uhL11', 'uhR11')):
        test = i.express(AhX).subs(vel_m_dict).simplify()
        print(j +  ' only components in AhX.x direction? ', sm.Abs(me.dot(test, AhX.y)) + sm.Abs(me.dot(test, AhX.z)) == 0)
    print('\n')
    '''
    # get this component
    uvL1 = me.dot(uvL11, AvX.x)
    uvR1 = me.dot(uvR11, AvX.x)
    uhL1 = me.dot(uhL11, AhX.x)
    uhR1 = me.dot(uhR11, AhX.x)
    
    # collect the speed constraints in a list for later use
    subs_dict_l = {sm.Derivative(q3, t): u3, sm.Derivative(q4, t): u4, sm.Derivative(q5, t): u5}
    loesung = [i.subs(subs_dict_l) for i in (vel_m_dict[u1], vel_m_dict[u2], uvL1, uvR1, uhL1, uhR1)]
    print('loesung DS', set().union(*[me.find_dynamicsymbols(loesung[i]) for i in range(len(loesung))]), '\n')
    
    # as it should be, these speeds only have components in AvX.z / AhX.z direction 
    print('speed(DmcvL) = ', ((DmcvL.vel(N).express(AvX).subs(vel_m_dict)).subs(subs_dict_l)).simplify())
    print('speed(DmchR) = ', ((DmchR.vel(N).express(AhX).subs(vel_m_dict)).subs(subs_dict_l)).simplify())

Define the various bodies.

As the energy of a system may give hints that maybe a mistake was made
setting up Kane’s equations, e.g. absent any friction or applied
forces/torques the total energy must be constant (this has helped me
MANY times!), I like to calculate them and look at them.

.. jupyter-execute::

    IvL = me.inertia(AvL, iXX, iYY, iZZ)
    IvR = me.inertia(AvR, iXX, iYY, iZZ)
    IhL = me.inertia(AhL, iXX, iYY, iZZ)
    IhR = me.inertia(AhR, iXX, iYY, iZZ)
    
    IvX = me.inertia(AvX, iXXa, iYYa, iZZa)
    IhX = me.inertia(AhX, iXXa, iYYa, iZZa)
    IV = me.inertia(AV, iXXV, iYYV, iZZV)
    
    BodyvL = me.RigidBody('BodyvL', DmcvL, AvL, m, (IvL, DmcvL))
    BodyvR = me.RigidBody('BodyvR', DmcvR, AvR, m, (IvR, DmcvR))
    BodyhL = me.RigidBody('BodyhL', DmchL, AhL, m, (IhL, DmchL))
    BodyhR = me.RigidBody('BodyhR', DmchR, AhR, m, (IhR, DmchR))
    
    BodyvX = me.RigidBody('BodyvX', mv, AvX, ma, (IvX, mv))
    BodyhX = me.RigidBody('BodyhX', mh, AhX, ma, (IhX, mh))
    BodyV = me.RigidBody('BodyV', mm, AV, mV, (IV, mm))
    
    observervL = me.Particle('observervL', ovL, mo)
    observervR = me.Particle('observervR', ovR, mo)
    observerhL = me.Particle('observerhL', ohL, mo)
    observerhR = me.Particle('observerhR', ohR, mo)
    
    BODY = [BodyvL, BodyvR, BodyhL, BodyhR, BodyvX, BodyhX, BodyV, observervL, observervR, observerhL, observerhR]
    
    punkte = [DmcvL, DmcvR, DmchL, DmchR, mv, mm, mh, ovL, ovR, ohL, ohR]
    massen = [m] * 4 + [ma, mV, ma]   + [mo] * 4
    
    pot_energie = sum([i * g * me.dot(j.pos_from(P0), N.y) for i, j in zip(massen, punkte)])
    
    subs_dict_kin = {sm.Derivative(i, t): j for i, j in zip((q1, q2, q3, q4, q5, qvL, qvR, qhL, qhR), (u1, u2, u3, u4, u5, uvL, uvR, uhL, uhR))}
    kin_energie = (sum([koerper.kinetic_energy(N) for koerper in BODY]).subs(subs_dict_kin))
    
    print('kinetic energy DS:', me.find_dynamicsymbols(kin_energie))
    print('kinetic energy free symbols', kin_energie.free_symbols)
    print('pot. energy DS:', me.find_dynamicsymbols(pot_energie))

Set up the external forces acting on the system

.. jupyter-execute::

    FL1 = [(i, -j*g*N.y) for i, j in zip(punkte, massen)]
    raeder = [AvL, AvR, AhL, AhR]
    raeder_speed = [uvL, uvR, uhL, uhR]
    Torque = [(i, -reibung * j * i.x) for i, j in zip(raeder, raeder_speed)]
    FL = FL1 + Torque

Set up **Kane’s equations** of motion. Following JM’s example, I use
:math:`u_3, u_4, u_5` as independent speeds,
:math:`u_1, u_2, u_{vL}, u_{vR}, u_{uhL}, u_{uhR}` as dependent speeds.

I find it useful, do look at the dynamic symbols /free symbols of the
mass matrix and of the force. This often helps to correct errors in
setting up Kane’s equations. Looking at the number of operations gives a
feel how ‘large’ the problem is.

.. jupyter-execute::

    q_ind = [q1, q2, q3, q4, q5] + [qvL, qvR, qhL, qhR]
    u_ind = [u3, u4, u5]
    u_dep = [u1, u2, uvL, uvR, uhL, uhR] 
    
    kd = [i - j.diff(t) for i, j in zip((u1, u2, u3, u4, u5, uvL, uvR, uhL, uhR), q_ind)]
    speed_constraint = [i - j for i, j in zip(loesung, u_dep)]
    
    KM = me.KanesMethod(N, q_ind=q_ind, u_ind=u_ind, u_dependent=u_dep, kd_eqs=kd, velocity_constraints=speed_constraint)
    (fr, frstar) = KM.kanes_equations(BODY, FL)
    
    MM = KM.mass_matrix_full
    print('MM DS', me.find_dynamicsymbols(MM))
    print('MM free symbols', MM.free_symbols)
    print('MM has {} operations'.format(sum([MM[i, j].count_ops(visual=False) 
        for i in range(MM.shape[0]) for j in range(MM.shape[1])])), '\n')
    
    force = KM.forcing_full
    print('shape of force', force.shape)
    print('force DS', me.find_dynamicsymbols(force))
    print('force free symbols', force.free_symbols)
    print('force has {} operations'.format(sum([force[i].count_ops(visual=False) 
                    for i in range(len(force))])), '\n')

All these points / locations of points are needed later for the
animation only.

.. jupyter-execute::

    Dmc_pos = [me.dot(i.pos_from(P0), uv) for i in (DmcvL, DmcvR, DmchL, DmchR) for uv in (N.x, N.z)]
    observer_pos = [me.dot(i.pos_from(P0), uv) for i in (ovL, ovR, ohL, ohR) for uv in (N.x, N.z)]
    mitte_pos = [me.dot(i.pos_from(P0), uv) for i in (mv, mm, mh) for uv in (N.x, N.z)]
    
    radvL1, radvL2, radvR1, radvR2 = sm.symbols('radvL1, radvL2, radvR1, radvR2', cls=me.Point) # just for the wheels in animation
    radhL1, radhL2, radhR1, radhR2 = sm.symbols('radhL1, radhL2, radhR1, radhR2', cls=me.Point) # just for the wheels in animation
    
    radvL1.set_pos(DmcvL, -r*AvX.z)
    radvL2.set_pos(DmcvL, r*AvX.z)
    radvR1.set_pos(DmcvR, -r*AvX.z)
    radvR2.set_pos(DmcvR, r*AvX.z)
    
    radhL1.set_pos(DmchL, -r*AhX.z)
    radhL2.set_pos(DmchL, r*AhX.z)
    radhR1.set_pos(DmchR, -r*AhX.z)
    radhR2.set_pos(DmchR, r*AhX.z)
    
    rad_liste = [radvL1, radvL2, radvR1, radvR2, radhL1, radhL2, radhR1, radhR2]
    rad_pos = [me.dot(i.pos_from(P0), uv) for i in rad_liste for uv in (N.x, N.z)]

Convert all the sympy functions to numpy functions, so numerical
calculations may be done. *cse=True* speeds up the numerical integration
substantially.

.. jupyter-execute::

    pL = [m, mo, ma, mV, g, r, l1, l2,  iXX, iYY, iZZ, iXXa, iYYa, iZZa, iXXV, iYYV, iZZV, reibung, alphavL, betavL, alphavR, betavR, alphahL, betahL, alphahR, betahR ]
    print('number of parameters:', len(pL))
    qL1 = q_ind + u_ind + u_dep
    print('vector of variables is:', qL1, ', its length is', len(qL1), '\n')
    MM_lam = sm.lambdify(qL1 + pL, MM, cse=True)
    force_lam = sm.lambdify(qL1 + pL, force, cse=True)
    
    pot_lam = sm.lambdify(qL1 + pL, pot_energie, cse=True)
    kin_lam = sm.lambdify(qL1 + pL, kin_energie, cse=True)
    
    Dmc_lam = sm.lambdify(qL1 + pL, Dmc_pos, cse=True)
    observer_lam = sm.lambdify(qL1 + pL, observer_pos, cse=True)
    mitte_lam = sm.lambdify(qL1 + pL, mitte_pos, cse=True)
    rad_lam = sm.lambdify(qL1 + pL, rad_pos, cse=True)
    
    loesung_lam = sm.lambdify([q1, q2, q3, q4, q5, u3, u4, u5] + pL, loesung, cse=True)
    
    print('it took {:.3f} sec to establish Kanes equations'.format(time.time() - start))

Perform the **numerical integration**. While it makes sense to name the
input variables similar to their sympy symbols / dynamic symbols
equivalents, *avoid using the*\ **same**\ *names*. This would overwrite
these symbols, with unpredictable consequences.

Input data: - :math:`m_1`: mass of the ball - :math:`m_{o1}`: mass of
the particle (observer) - :math:`m_{a1}` mass of axles - :math:`m_{V1}`
mass of connecting axle -
:math:`\alpha_{L1}, \beta_{L1}, \alpha_{R1}, \beta_{R1}`: define the
location of the particle relative to the center of the ball, 0 <=
:math:`\alpha_{...} , \beta_{...}` < 1. - :math:`q_{31},..., q_{51}`:
initial generalized coordinates - :math:`u_{21},..., u_{51}`: initial
speeds - intervall: the integration will run from 0 to intervall -
schritte: number of time instances returned by ivp_solve

Making :math:`m_o` very small results (obviously) in a motion very
similar to the one in JM’s lecture. Making in ‘similar’ to :math:`m`
results in ‘wilder’ motions.

.. jupyter-execute::

    #==============================================================
    start2 = time.time()
    # set the input variables.
    m1 = 1. 
    mo1 = 5.e-2
    ma1 = 1.
    mV1 = 10.
    l11 = 2.
    l21 = 4.
    r1 = 1.
    reibung1 = 0.
    
    alphavL1, betavL1, alphavR1, betavR1 = 0.999, 0., 0., 0.999
    alphahL1, betahL1, alphahR1, betahR1 = 0., 0.999, 0., 0.999
    
    q11, q21, q31, q41, q51 = 0., 0., np.pi/2. + 0.1, -np.deg2rad(5), np.deg2rad(5)
    qvL1, qvR1, qhL1, qhR1 = 0., 0., 0., 0.
    u31, u41, u51 = 0.05, 0.005, -0.005
    
    intervall = 50.
    #================================================================
    schritte = int(intervall * 30)
    print('Arguments:')
    print('[m, mo, ma, mV, g, r, l1, l2,  iXX, iYY, iZZ, iXXa, iYYa, iZZa, iXXV, iYYV, iZZV, reibung, alphavL, betavL, alphavR, betavR, alphahL, betahL, alphahR, betahR ] = ')
    iXX1 = 0.5 * m1 * r1**2   # from the internet
    iYY1 = 0.25 * m1 * r1**2
    iZZ1 = iYY1
    
    iXXa1 = 0.
    iYYa1 = 1./12. * ma1 * l11**2   # from the internet for a thin rod
    iZZa1 =0.
    
    iXXV1 = 0.
    iYYV1 = 1./12. * mV1 * l21**2
    iZZV1 =0.
    
    pL_vals = [m1, mo1, ma1, mV1, 9.81, r1, l11, l21, iXX1, iYY1, iZZ1, iXXa1, iYYa1, iZZa1, iXXV1, iYYV1, iZZV1, reibung1, alphavL1, betavL1, alphavR1, betavR1, alphahL1, betahL1, alphahR1, betahR1]
    print(pL_vals)
    print('number of parameters = ', len(pL_vals))

Find the initial values of the “dependent” variables: -
:math:`u_1, u_2, u_{vL}, u_{vR}, u_{hL} u_{hR}`

As the wheel is rotationally symmetric the dependent variables
:math:`q_{vL}, q_{vR}, q_{hL}, q_{hR}` may have any initial values. I
arbitrarily set them to zero above.

.. jupyter-execute::

    u11, u21, uvL111, uvR111, uhL111, uhR111 = loesung_lam(q11, q21, q31, q41, q51, u31, u41, u51, *pL_vals)
    print('u11, u21, uvL1, uvR1, uhL1, uhR1 = ', ['{:.3f}'.format(i) for i in (u11, u21, uvL111, uvR111, uhL111, uhR111)])
    
    y0 = [q11, q21, q31, q41, q51, qvL1, qvR1, qhL1, qhR1] + [u31, u41, u51, u11, u21, uvL111, uvR111, uhL111, uhR111]
    
    print('starting vector: ', ['{:.3f}'.format(y0[i]) for i in range(len(y0))])

Ensure that the particle (observer) is inside the ball, that is $
:raw-latex:`\alpha`^2 + :raw-latex:`\beta`^2 + :raw-latex:`\gamma`^2 <
1. $ If this is not the case, an Exception is raised.

.. jupyter-execute::

    if alphavL1**2 + betavL1**2 >= 1. or alphavR1**2 + betavR1**2 >= 1. or alphahL1**2 + betahL1**2 >=1 or alphahR1**2 + betahR1**2 >= 1.:
        raise Exception('the observer on at least one wheel is outside of the wheel')

**Numerical integration**.

You may set atol, rtol higher than the standard, to get better results,
that is total energy closer to being constant absent any friction. Of
course at the expense of longer duration of the calculation. *method =
‘Radau’* seems to give better constant total energy in absence of
friction than no method. solve_ivp gives messages as to how the
integrationworked, I print them.

.. jupyter-execute::

    times = np.linspace(0, intervall, schritte)
      
    def gradient(t, y, args):
        vals = np.concatenate((y, args))
        sol = np.linalg.solve(MM_lam(*vals), force_lam(*vals))
        return np.array(sol).T[0]
    
    t_span = (0., intervall)
    resultat1 = solve_ivp(gradient, t_span, y0, t_eval=times, args=(pL_vals,), method='Radau') #, atol=1.e-6, rtol=1.e-6)
    resultat = resultat1.y.T
    print(resultat.shape)
    event_dict = {-1: 'Integration failed', 0: 'Integration finished successfully', 1: 'some termination event'}
    print(event_dict[resultat1.status], ', message was:', resultat1.message)
        
    print('resultat shape', resultat.shape, '\n')
    
    print("To numerically integrate an intervall of {} sec, the routine cycled {} times and it took {:.3f} sec "
          .format(intervall, resultat1.nfev, time.time() - start2))

Plot whichever coordinates of the system you may want to see.

.. jupyter-execute::

    bezeichnung = ['q11', 'q21', 'q31', 'q41', 'q51', 'qvL1', 'qvR1','qhL1', 'qhR1'] +  ['u11', 'u21', 'u31', 'u41', 'u51', 'uvL1', 'uvR1', 'uhL1', 'uhR1']
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in (14, 15, 16, 17):
        ax.plot(times[: resultat.shape[0]], resultat[:, i], label = bezeichnung[i])
    ax.set_title('generalized coordinates')
    ax.set_xlabel('time (sec)')
    ax.legend();

Plot the energies of the ball. Absent any friction, the total energy
must be constant.

.. jupyter-execute::

    pot_np = np.empty(resultat.shape[0])
    kin_np = np.empty(resultat.shape[0])
    total_np = np.empty(resultat.shape[0])
    for l in range(resultat.shape[0]):
        pot_np[l] = pot_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)
        kin_np[l] = kin_lam(*[resultat[l, j] for j in range(resultat.shape[1])], *pL_vals)
        total_np[l] = pot_np[l] + kin_np[l]
    
    if reibung1 == 0.:
        fehler = ((x:=max(total_np)) - min(total_np)) / x * 100.
        print('max. deviation from total energy = constant is {:.3f} % of max total energy'.format(fehler))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times[: resultat.shape[0]], kin_np, label = 'kin energy')
    ax.plot(times[: resultat.shape[0]], pot_np, label = 'pot energy')
    ax.plot(times[: resultat.shape[0]], total_np, label = 'total energy')
    ax.set_title('Energy of the system, friction = {}'.format(reibung1))
    ax.set_xlabel('time (sec)')
    ax.set_ylabel('energy (Nm)')
    ax.legend();

**Animation**

Movement of the snakeboard in the X / Z direction. The small dots on the
wheels show the projection into the X / Z plane of the locations of the
particles mounted to the wheels. HTML is needed to show the animation on
my iPad. no idea whether needed on other machines. It is SLOW!

.. jupyter-execute::

    # to enable animtion to run, even in the integration did not finish
    schritte1 = schritte
    schritte = resultat.shape[0]
    times1 = times
    times = times1[: schritte]
    
    Dmc_list = []
    observer_list = []
    mitte_list = []
    rad_list = []
    
    for i in range(schritte):
        Dmc_list.append(Dmc_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals))
        observer_list.append(observer_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals))
        mitte_list.append(mitte_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals))
        rad_list.append(rad_lam(*[resultat[i, j] for j in range(resultat.shape[1])], *pL_vals))
    
    # Determine the size of the picture
    xmin = min([Dmc_list[i][j] for i in range(schritte) for j in range(8)])
    xmax = max([Dmc_list[i][j] for i in range(schritte) for j in range(8)])
    
    def animate_pendulum(times, Dmc_list, mitte_list, observer_list, rad_list):    
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'aspect': 'equal'})
        ax.axis('on')
        ax.set_xlim(xmin - 1., xmax + 1.)
        ax.set_ylim(xmin - 1., xmax + 1.)
        ax.set_xlabel(' X Axis', fontsize=15)
        ax.set_ylabel('Z Axis', fontsize=15)
    
        line1, = ax.plot([], [], 'o', markersize=10, color='green') # DmcvL     
        line2, = ax.plot([], [], 'o', markersize=10, color='red')   # DmcvR
        line3, = ax.plot([], [], 'o', markersize=10, color='black') # DmchL
        line4, = ax.plot([], [], 'o', markersize=10, color='violet') # DmchR
     
        line5, = ax.plot([], [], color='green', linewidth=0.25)        # tracing the motion of DmcvL
        line6, = ax.plot([], [], color='red', linewidth=0.25)         # dto for DmcvR
        line7, = ax.plot([], [], color='black', linewidth=0.25)        # tracing the motion of DmhL
        line8, = ax.plot([], [], color='violet', linewidth=0.25)         # dto for DmhcR
    
        line9, = ax.plot([], [], 'bo', linestyle='-', linewidth=2., markersize=0.)  # the axle front
        line10, = ax.plot([], [], 'bo', linestyle='-', linewidth=2., markersize=0.)  # the axle rear
        line11, = ax.plot([], [], 'bo', linestyle='-', linewidth=2., markersize=0.)  # the connecting axle
        
        line12, = ax.plot([], [], 'bo', linestyle='-', linewidth=1., markersize=0.)  # front left wheel
        line13, = ax.plot([], [], 'bo', linestyle='-', linewidth=1., markersize=0.)  
        line14, = ax.plot([], [], 'bo', linestyle='-', linewidth=1., markersize=0.)  
        line15, = ax.plot([], [], 'bo', linestyle='-', linewidth=1., markersize=0.)  
        
        line16, = ax.plot([], [], 'o', markersize=5, color='green') # observervL     
        line17, = ax.plot([], [], 'o', markersize=5, color='red')   
        line18, = ax.plot([], [], 'o', markersize=5, color='black')
        line19, = ax.plot([], [], 'o', markersize=5, color='violet')
        
        def animate(i):
            ax.set_title('Running time {:.1f} sec'.format(times[i]), fontsize=15)
            line1.set_data(Dmc_list[i][0], Dmc_list[i][1])      # DmcvL, etc
            line2.set_data(Dmc_list[i][2], Dmc_list[i][3])   
            line3.set_data(Dmc_list[i][4], Dmc_list[i][5])  
            line4.set_data(Dmc_list[i][6], Dmc_list[i][7])    
    
            line5.set_data([Dmc_list[j][0] for j in range(i+1) ], [Dmc_list[j][1] for j in range(i+1) ])      # tracing DmcvL, etc  
            line6.set_data([Dmc_list[j][2] for j in range(i+1) ], [Dmc_list[j][3] for j in range(i+1) ])      # tracing DmcvL, etc        
            line7.set_data([Dmc_list[j][4] for j in range(i+1) ], [Dmc_list[j][5] for j in range(i+1) ])      # tracing DmcvL, etc        
            line8.set_data([Dmc_list[j][6] for j in range(i+1) ], [Dmc_list[j][7] for j in range(i+1) ])      # tracing DmcvL, etc         
                
            x_values = [Dmc_list[i][0], Dmc_list[i][2]]           # axle front
            z_values = [Dmc_list[i][1], Dmc_list[i][3]]
            line9.set_data(x_values, z_values)
            
            x_values = [Dmc_list[i][4], Dmc_list[i][6]]           # axle rear
            z_values = [Dmc_list[i][5], Dmc_list[i][7]]
            line10.set_data(x_values, z_values)
            
            x_values = [mitte_list[i][0], mitte_list[i][4]]          # connecting axle
            z_values = [mitte_list[i][1], mitte_list[i][5]]
            line11.set_data(x_values, z_values)
            
            x_values = [rad_list[i][0], rad_list[i][2]]          # front left wheel
            z_values = [rad_list[i][1], rad_list[i][3]]
            line12.set_data(x_values, z_values)
            
            x_values = [rad_list[i][4], rad_list[i][6]]           
            z_values = [rad_list[i][5], rad_list[i][7]]
            line13.set_data(x_values, z_values)
            
            x_values = [rad_list[i][8], rad_list[i][10]]           
            z_values = [rad_list[i][9], rad_list[i][11]]
            line14.set_data(x_values, z_values)
            
            x_values = [rad_list[i][12], rad_list[i][14]]           
            z_values = [rad_list[i][13], rad_list[i][15]]
            line15.set_data(x_values, z_values)
            
            line16.set_data(observer_list[i][0], observer_list[i][1])      # DmcvL, etc
            line17.set_data(observer_list[i][2], observer_list[i][3])   
            line18.set_data(observer_list[i][4], observer_list[i][5])  
            line19.set_data(observer_list[i][6], observer_list[i][7])
                    
            return line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, 
        line12, line13, line14, line15, line16, line17, line18, line19,
    
        anim = animation.FuncAnimation(fig, animate, frames=schritte,
                                       interval=1000*max(times) / schritte,
                                       blit=True)
        plt.close(fig)
        return anim
    
    anim = animate_pendulum(times, Dmc_list, mitte_list, observer_list, rad_list)
    print('it took, before HTML started, {:.3f} sec to run the program'.format(time.time() - start))
    HTML(anim.to_jshtml())

