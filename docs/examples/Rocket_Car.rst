##

.. raw:: html

   <center>

:math:`\mathbf{\text{Modeling of a Variable-Mass Nonholonomic Gyrostatic Rocket Car Using Extended Kane's Equations}}`

.. raw:: html

   </center>

#####

.. raw:: html

   <center>

:math:`\text{Ge, Z.M. and Cheng, Y.H., 1982. Extended Kane’s equations for nonholonomic variable mass system.}`

.. raw:: html

   </center>

#####

.. raw:: html

   <center>

:math:`\text{Figure 1: Idealized model of a jet racing car that is propelled by a rocket engine at point P, the rocket engine being treated as a variable mass particle at P}`

.. raw:: html

   </center>

#####

.. raw:: html

   <center>

:math:`\text{Here}, \{a_x: g_3, a_y: g_1, a_z: g_2\}`

.. raw:: html

   </center>

#####

.. raw:: html

   <center>

:math:`\text{Leimanis, E., 2013. The general problem of the motion of coupled rigid bodies about a fixed point (Vol. 7). Springer Science & Business Media. (Page 207)}`

.. raw:: html

   </center>

.. code:: ipython3

    import sympy as sm
    import sympy.physics.mechanics as me
    from pydy.system import System
    import numpy as np
    from sympy.simplify.fu import TR2
    import matplotlib.pyplot as plt
    from scipy.integrate import odeint
    me.init_vprinting()

.. code:: ipython3

    N = me.ReferenceFrame('N')
    
    q1, q2, q3, q4, q5, q6, q7, q8 = me.dynamicsymbols('q1:9')
    
    A2 = N.orientnew('A_2', 'Axis', [q3, N.z])
    A1 = A2.orientnew('A_1', 'Axis', [q4, A2.z])
    
    B1 = A1.orientnew('B_1', 'Axis', [q5, A1.y])
    B2 = A1.orientnew('B_2', 'Axis', [q6, A1.y])
    B3 = A2.orientnew('B_3', 'Axis', [q7, A2.y])
    B4 = A2.orientnew('B_4', 'Axis', [q8, A2.y])
    
    t = me.dynamicsymbols._t

.. code:: ipython3

    O = me.Point('O') # fixed point in the inertial reference frame
    O.set_vel(N, 0)

.. code:: ipython3

    L, l , a, b, r1, r2 = sm.symbols('L, l , a, b, r_1, r_2')

.. code:: ipython3

    Q = O.locatenew('Q', q1 * N.x + q2 * N.y)

.. code:: ipython3

    P = Q.locatenew('P', L * -A2.x)

.. code:: ipython3

    C = P.locatenew('C', l * A2.x)

.. code:: ipython3

    Q.set_vel(N, Q.pos_from(O).dt(N))
    Q.vel(N)




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y}



.. code:: ipython3

    P.v2pt_theory(Q, N, A2)
    P.vel(N)




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} -  L \dot{q}_{3}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    C.v2pt_theory(P, N, A2)
    # C.vel(N)




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + (- L \dot{q}_{3} + l \dot{q}_{3})\mathbf{\hat{a_2}_y}



.. code:: ipython3

    A1.ang_vel_in(A2).express(A1)




.. math::

    \displaystyle \dot{q}_{4}\mathbf{\hat{a_1}_z}



.. code:: ipython3

    u1, u2 = me.dynamicsymbols('u_1:3')

.. code:: ipython3

    z1 = sm.Eq(u1, A1.ang_vel_in(A2).dot(A1.z))
    z2 = sm.Eq(u2, Q.vel(N).dot(A1.x))

.. code:: ipython3

    u = sm.trigsimp(sm.solve([z1, z2], A1.ang_vel_in(A2).dot(A1.z), Q.vel(N).dot(A1.x)))
    u




.. math::

    \displaystyle \left\{ \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{1} : u_{2}, \  \dot{q}_{4} : u_{1}\right\}



:math:`\text{Nonholonomic Constraints:}\ B_1`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    B1_center = Q.locatenew('B_1_center', a * A1.y)
    B1_center.pos_from(Q)




.. math::

    \displaystyle a\mathbf{\hat{a_1}_y}



.. code:: ipython3

    B1_center.v2pt_theory(Q, N, A1)
    B1_center.vel(N).express(A1).simplify()




.. math::

    \displaystyle (- a \left(\dot{q}_{3} + \dot{q}_{4}\right) + \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{1})\mathbf{\hat{a_1}_x} + (- \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{2})\mathbf{\hat{a_1}_y}



.. code:: ipython3

    B1_ground = B1_center.locatenew('B_1_ground', r1 * -A1.z)
    B1_ground.pos_from(B1_center)




.. math::

    \displaystyle -  r_{1}\mathbf{\hat{a_1}_z}



.. code:: ipython3

    B1_ground.v2pt_theory(B1_center, N, B1)
    B1_ground.vel(N).simplify()




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + (- a \left(\dot{q}_{3} + \dot{q}_{4}\right) - r_{1} \dot{q}_{5})\mathbf{\hat{a_1}_x}



.. code:: ipython3

    B1_cons = [me.dot(B1_ground.vel(N).simplify(), uv) for uv in A1]
    sm.trigsimp(B1_cons)




.. math::

    \displaystyle \left[ - a \left(\dot{q}_{3} + \dot{q}_{4}\right) - r_{1} \dot{q}_{5} + \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right)\right) \dot{q}_{1} + \left(\operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{q}_{2}, \  \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right)\right) \dot{q}_{2} + \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right) - \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{q}_{1}, \  0\right]



.. code:: ipython3

    eq1 = sm.Eq(B1_cons[0].simplify().subs(u), 0)
    eq1




.. math::

    \displaystyle - a \left(u_{1} + \dot{q}_{3}\right) - r_{1} \dot{q}_{5} + u_{2} = 0



.. code:: ipython3

    eq2 = sm.Eq(B1_cons[1].simplify().subs(u), 0)
    eq2




.. math::

    \displaystyle - \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{2} = 0



:math:`\text{Nonholonomic Constraints:}\ B_2`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    B2_center = Q.locatenew('B_1_center', a * -A1.y)
    B2_center.pos_from(Q)




.. math::

    \displaystyle -  a\mathbf{\hat{a_1}_y}



.. code:: ipython3

    B2_center.v2pt_theory(Q, N, A1)
    B2_center.vel(N).express(A1).simplify()




.. math::

    \displaystyle (a \left(\dot{q}_{3} + \dot{q}_{4}\right) + \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{1})\mathbf{\hat{a_1}_x} + (- \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{2})\mathbf{\hat{a_1}_y}



.. code:: ipython3

    B2_ground = B2_center.locatenew('B_2_ground', r1 * -A1.z)
    B2_ground.pos_from(B2_center)




.. math::

    \displaystyle -  r_{1}\mathbf{\hat{a_1}_z}



.. code:: ipython3

    B2_ground.v2pt_theory(B2_center, N, B2)
    B2_ground.vel(N).simplify()




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + (a \left(\dot{q}_{3} + \dot{q}_{4}\right) - r_{1} \dot{q}_{6})\mathbf{\hat{a_1}_x}



.. code:: ipython3

    B2_cons = [me.dot(B2_ground.vel(N).simplify(), uv) for uv in A1]
    sm.trigsimp(B2_cons)




.. math::

    \displaystyle \left[ a \left(\dot{q}_{3} + \dot{q}_{4}\right) - r_{1} \dot{q}_{6} + \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right)\right) \dot{q}_{1} + \left(\operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{q}_{2}, \  \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right)\right) \dot{q}_{2} + \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right) - \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{q}_{1}, \  0\right]



.. code:: ipython3

    eq3 = sm.Eq(B2_cons[0].simplify().subs(u), 0)
    eq3




.. math::

    \displaystyle a \left(u_{1} + \dot{q}_{3}\right) - r_{1} \dot{q}_{6} + u_{2} = 0



.. code:: ipython3

    eq4 = sm.Eq(B2_cons[1].simplify().subs(u), 0)
    eq4




.. math::

    \displaystyle - \operatorname{sin}\left(q_{3} + q_{4}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3} + q_{4}\right) \dot{q}_{2} = 0



:math:`\text{Nonholonomic Constraints:}\ B_3`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    B3_center = P.locatenew('B_3_center', b * A2.y)
    B3_center.pos_from(P)




.. math::

    \displaystyle b\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B3_center.v2pt_theory(P, N, A2)
    B3_center.vel(N).express(A2).simplify()




.. math::

    \displaystyle (- b \dot{q}_{3} + \operatorname{sin}\left(q_{3}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{1})\mathbf{\hat{a_2}_x} + (- L \dot{q}_{3} - \operatorname{sin}\left(q_{3}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{2})\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B3_ground = B3_center.locatenew('B_3_ground', r2 * -A2.z)
    B3_ground.pos_from(B3_center)




.. math::

    \displaystyle -  r_{2}\mathbf{\hat{a_2}_z}



.. code:: ipython3

    B3_ground.v2pt_theory(B3_center, N, B3)
    B3_ground.vel(N).simplify()




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + (- b \dot{q}_{3} - r_{2} \dot{q}_{7})\mathbf{\hat{a_2}_x} -  L \dot{q}_{3}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B3_cons = [me.dot(B3_ground.vel(N).simplify(), uv) for uv in A2]
    sm.trigsimp(B3_cons)




.. math::

    \displaystyle \left[ - b \dot{q}_{3} - r_{2} \dot{q}_{7} + \operatorname{sin}\left(q_{3}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{1}, \  - L \dot{q}_{3} - \operatorname{sin}\left(q_{3}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{2}, \  0\right]



.. code:: ipython3

    eq5 = sm.Eq(B3_cons[0].simplify().subs(u), 0)
    eq5




.. math::

    \displaystyle - b \dot{q}_{3} - r_{2} \dot{q}_{7} + \operatorname{sin}\left(q_{3}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{1} = 0



.. code:: ipython3

    eq6 = sm.Eq(B3_cons[1].simplify().subs(u), 0)
    eq6




.. math::

    \displaystyle - L \dot{q}_{3} - \operatorname{sin}\left(q_{3}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{2} = 0



:math:`\text{Nonholonomic Constraints:}\ B_4`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    B4_center = P.locatenew('B_4_center', b * -A2.y)
    B4_center.pos_from(P)




.. math::

    \displaystyle -  b\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B4_center.v2pt_theory(P, N, A2)
    # B4_center.vel(N).express(A2).simplify() # Invalid Json




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + b \dot{q}_{3}\mathbf{\hat{a_2}_x} -  L \dot{q}_{3}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B4_ground = B4_center.locatenew('B_4_ground', r2 * -A2.z)
    B4_ground.pos_from(B4_center)




.. math::

    \displaystyle -  r_{2}\mathbf{\hat{a_2}_z}



.. code:: ipython3

    B4_ground.v2pt_theory(B4_center, N, B4)
    B4_ground.vel(N).simplify()




.. math::

    \displaystyle \dot{q}_{1}\mathbf{\hat{n}_x} + \dot{q}_{2}\mathbf{\hat{n}_y} + (b \dot{q}_{3} - r_{2} \dot{q}_{8})\mathbf{\hat{a_2}_x} -  L \dot{q}_{3}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    B4_cons = [me.dot(B4_ground.vel(N).simplify(), uv) for uv in A2]
    sm.trigsimp(B4_cons)




.. math::

    \displaystyle \left[ b \dot{q}_{3} - r_{2} \dot{q}_{8} + \operatorname{sin}\left(q_{3}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{1}, \  - L \dot{q}_{3} - \operatorname{sin}\left(q_{3}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{2}, \  0\right]



.. code:: ipython3

    eq7 = sm.Eq(B4_cons[0].simplify().subs(u), 0)
    eq7




.. math::

    \displaystyle b \dot{q}_{3} - r_{2} \dot{q}_{8} + \operatorname{sin}\left(q_{3}\right) \dot{q}_{2} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{1} = 0



.. code:: ipython3

    eq8 = sm.Eq(B4_cons[1].simplify().subs(u), 0)
    eq8




.. math::

    \displaystyle - L \dot{q}_{3} - \operatorname{sin}\left(q_{3}\right) \dot{q}_{1} + \operatorname{cos}\left(q_{3}\right) \dot{q}_{2} = 0



:math:`\text{LHS} \Longleftrightarrow \text{RHS}\ \text{in}\ z_1, z_2 \rightarrow \text{Equations}\ 9, 10`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    eq9 = sm.Eq(A1.ang_vel_in(A2).dot(A1.z), u1)
    eq9




.. math::

    \displaystyle \dot{q}_{4} = u_{1}



.. code:: ipython3

    eq10 = sm.Eq(Q.vel(N).dot(A1.x), u2)
    eq10




.. math::

    \displaystyle \left(- \operatorname{sin}\left(q_{3}\right) \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right)\right) \dot{q}_{1} + \left(\operatorname{sin}\left(q_{3}\right) \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{3}\right)\right) \dot{q}_{2} = u_{2}



:math:`\text{Solution of the System of Linear (in}\ \dot{q}_1,...\text{) Equations}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    solution = sm.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10], q1.diff(), q2.diff(), q3.diff(),  q4.diff(), q5.diff(), q6.diff(), q7.diff(), q8.diff())

.. code:: ipython3

    solution




.. math::

    \displaystyle \left\{ \dot{q}_{1} : u_{2} \operatorname{cos}\left(q_{3} + q_{4}\right), \  \dot{q}_{2} : u_{2} \operatorname{sin}\left(q_{3} + q_{4}\right), \  \dot{q}_{3} : \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  \dot{q}_{4} : u_{1}, \  \dot{q}_{5} : \frac{- L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) \left(a u_{1} + u_{2}\right) + 2 L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) u_{2} - a \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L r_{1} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  \dot{q}_{6} : \frac{L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) \left(a u_{1} + u_{2}\right) + a \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L r_{1} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  \dot{q}_{7} : \frac{\left(L - b \operatorname{tan}\left(q_{4}\right)\right) \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  \dot{q}_{8} : \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}\right\}



.. code:: ipython3

    sollist_keys = list(solution.keys())
    sollist_keys




.. math::

    \displaystyle \left[ \dot{q}_{8}, \  \dot{q}_{4}, \  \dot{q}_{7}, \  \dot{q}_{6}, \  \dot{q}_{1}, \  \dot{q}_{5}, \  \dot{q}_{2}, \  \dot{q}_{3}\right]



.. code:: ipython3

    sollist_values = list(solution.values())
    sollist_values




.. math::

    \displaystyle \left[ \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}, \  u_{1}, \  \frac{\left(L - b \operatorname{tan}\left(q_{4}\right)\right) \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  \frac{L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) \left(a u_{1} + u_{2}\right) + a \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L r_{1} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  u_{2} \operatorname{cos}\left(q_{3} + q_{4}\right), \  \frac{- L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) \left(a u_{1} + u_{2}\right) + 2 L \left(L + b \operatorname{tan}\left(q_{4}\right)\right) u_{2} - a \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L r_{1} \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}, \  u_{2} \operatorname{sin}\left(q_{3} + q_{4}\right), \  \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2} \operatorname{tan}\left(q_{4}\right)}{L \left(L + b \operatorname{tan}\left(q_{4}\right)\right)}\right]



.. code:: ipython3

    sollist_values_simple = []
    for i in range(8):
       sollist_values_simple.append(sm.factor(TR2(sollist_values[i]).simplify()))

.. code:: ipython3

    sollist_values_simple




.. math::

    \displaystyle \left[ \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}, \  u_{1}, \  - \frac{\left(- L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}, \  \frac{L a u_{1} + L u_{2} + a u_{2} \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  u_{2} \operatorname{cos}\left(q_{3} + q_{4}\right), \  - \frac{L a u_{1} - L u_{2} + a u_{2} \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  u_{2} \operatorname{sin}\left(q_{3} + q_{4}\right), \  \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L}\right]



.. code:: ipython3

    soldict = dict(zip(sollist_keys, sollist_values_simple)) 
    soldict




.. math::

    \displaystyle \left\{ \dot{q}_{1} : u_{2} \operatorname{cos}\left(q_{3} + q_{4}\right), \  \dot{q}_{2} : u_{2} \operatorname{sin}\left(q_{3} + q_{4}\right), \  \dot{q}_{3} : \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L}, \  \dot{q}_{4} : u_{1}, \  \dot{q}_{5} : - \frac{L a u_{1} - L u_{2} + a u_{2} \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  \dot{q}_{6} : \frac{L a u_{1} + L u_{2} + a u_{2} \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  \dot{q}_{7} : - \frac{\left(- L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}, \  \dot{q}_{8} : \frac{\left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) u_{2}}{L r_{2}}\right\}



:math:`\text{Reformulated Velocity and Angular Velocity Expressions}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    N_v_Q = Q.vel(N).subs(soldict).express(A1).simplify()
    N_v_Q




.. math::

    \displaystyle u_{2}\mathbf{\hat{a_1}_x}



.. code:: ipython3

    N_v_P = P.vel(N).subs(soldict).express(A2).simplify()
    N_v_P




.. math::

    \displaystyle u_{2} \operatorname{cos}\left(q_{4}\right)\mathbf{\hat{a_2}_x}



.. code:: ipython3

    N_v_C = C.vel(N).subs(soldict).express(A2).simplify()
    N_v_C




.. math::

    \displaystyle u_{2} \operatorname{cos}\left(q_{4}\right)\mathbf{\hat{a_2}_x} + \frac{l u_{2} \operatorname{sin}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    N_w_A1 = A1.ang_vel_in(N).subs(soldict).express(A1).simplify()
    N_w_A1




.. math::

    \displaystyle (u_{1} + \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L})\mathbf{\hat{a_1}_z}



.. code:: ipython3

    N_w_A2 = A2.ang_vel_in(N).subs(soldict).express(A2).simplify()
    N_w_A2




.. math::

    \displaystyle \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_z}



:math:`\text{Partial Velocities and Partial Angular Velocities}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    V_1_Q = N_v_Q.diff(u1, N)
    V_1_Q




.. math::

    \displaystyle 0



.. code:: ipython3

    V_2_Q = N_v_Q.diff(u2, N)
    V_2_Q




.. math::

    \displaystyle \mathbf{\hat{a_1}_x}



.. code:: ipython3

    V_1_C = N_v_C.diff(u1, N)
    V_1_C




.. math::

    \displaystyle 0



.. code:: ipython3

    V_2_C = N_v_C.diff(u2, N)
    V_2_C




.. math::

    \displaystyle \operatorname{cos}\left(q_{4}\right)\mathbf{\hat{a_2}_x} + \frac{l \operatorname{sin}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    V_1_P = N_v_P.diff(u1, N)
    V_1_P




.. math::

    \displaystyle 0



.. code:: ipython3

    V_2_P = N_v_P.diff(u2, N)
    V_2_P




.. math::

    \displaystyle \operatorname{cos}\left(q_{4}\right)\mathbf{\hat{a_2}_x}



.. code:: ipython3

    w_1_A1 = N_w_A1.diff(u1, N)
    w_1_A1




.. math::

    \displaystyle \mathbf{\hat{a_1}_z}



.. code:: ipython3

    w_2_A1 = N_w_A1.diff(u2, N)
    w_2_A1




.. math::

    \displaystyle \frac{\operatorname{sin}\left(q_{4}\right)}{L}\mathbf{\hat{a_1}_z}



.. code:: ipython3

    w_1_A2 = N_w_A2.diff(u1, N)
    w_1_A2




.. math::

    \displaystyle 0



.. code:: ipython3

    w_2_A2 = N_w_A2.diff(u2, N)
    w_2_A2




.. math::

    \displaystyle \frac{\operatorname{sin}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_z}



:math:`\text{Accelerations and Angular Accelerations}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    a_1__P, a_2__P, a_3__P, a_1__C, a_2__C, a_3__C, a__Q, alpha__A1, alpha__A2 = sm.symbols('a_1__P, a_2__P, a_3__P, a_1__C, a_2__C, a_3__C, a__Q, alpha__A1, alpha__A2')

.. code:: ipython3

    N_a_P = N_v_P.dt(N).subs(soldict)
    N_a_P




.. math::

    \displaystyle (- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2})\mathbf{\hat{a_2}_x} + \frac{u^{2}_{2} \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_y}



.. code:: ipython3

    N_a_C = N_v_C.dt(N).subs(soldict)
    N_a_C




.. math::

    \displaystyle (- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2} - \frac{l u^{2}_{2} \operatorname{sin}^{2}\left(q_{4}\right)}{L^{2}})\mathbf{\hat{a_2}_x} + (\frac{l u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{l \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L} + \frac{u^{2}_{2} \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{4}\right)}{L})\mathbf{\hat{a_2}_y}



.. code:: ipython3

    N_a_Q = N_v_Q.dt(N).subs(soldict)
    N_a_Q




.. math::

    \displaystyle \dot{u}_{2}\mathbf{\hat{a_1}_x} + \left(u_{1} + \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L}\right) u_{2}\mathbf{\hat{a_1}_y}



.. code:: ipython3

    N_aa_A1 = N_w_A1.dt(N).subs(soldict)
    N_aa_A1




.. math::

    \displaystyle (\dot{u}_{1} + \frac{u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{\operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L})\mathbf{\hat{a_1}_z}



.. code:: ipython3

    N_aa_A2 = N_w_A2.dt(N).subs(soldict)
    N_aa_A2




.. math::

    \displaystyle (\frac{u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{\operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L})\mathbf{\hat{a_2}_z}



:math:`\text{Forces}`
~~~~~~~~~~~~~~~~~~~~~

:math:`(F_r^*)_G = (F_r^*)_{GR} + (F_r^*)_{GI}`

where,

:math:`(F_r^*)_{GR} = {V_r}^G \cdot {F_G}^* + \omega_r^A \cdot {T_G}^*`

:math:`F_G^* = -m_G {a^G}^*`

:math:`T_G^* \overset{\Delta}{=} -[\alpha_A \cdot I_G + \omega_r^A \times (I_G \cdot \omega_r^A)]`

:math:`({F_r}^*)_{GI} = -J\{\omega_r^A [\ddot{q_k} g_1 + \dot{q_k} (\omega_3^A g_2 - \omega_2^A g_3)] + C_{kr} (\dot{\omega_1^A} + \ddot{q_k}) \}`

#####

.. raw:: html

   <center>

:math:`\text{Kane, T.R., 1978. Nonholonomic multibody systems containing gyrostats. In Dynamics of Multibody Systems (pp. 97-107). Springer, Berlin, Heidelberg.}`

.. raw:: html

   </center>

:math:`\text{Naming Convention:}`

:math:`({F_r}^*)_{G_n R}\ \text{(rigid)}`

:math:`({F_r}^*)_{G_n I}\ \text{(internal)}`

:math:`\text{Masses and Moments of Inertia}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    M1, M2 = sm.symbols('M_1, M_2')
    m = me.dynamicsymbols('m')

.. code:: ipython3

    I1x, I1y, I1z = sm.symbols('I_{1_x}, I_{1_y}, I_{1_z}')
    I2x, I2y, I2z = sm.symbols('I_{2_x}, I_{2_y}, I_{2_z}')
    J1, J2 = sm.symbols('J_1, J_2')

.. code:: ipython3

    I1 = me.inertia(A1, I1x, I1y, I1z)
    I1




.. math::

    \displaystyle I_{1_x}\mathbf{\hat{a_1}_x}\otimes \mathbf{\hat{a_1}_x} + I_{1_y}\mathbf{\hat{a_1}_y}\otimes \mathbf{\hat{a_1}_y} + I_{1_z}\mathbf{\hat{a_1}_z}\otimes \mathbf{\hat{a_1}_z}



.. code:: ipython3

    I2 = me.inertia(A2, I2x, I2y, I2z)
    I2




.. math::

    \displaystyle I_{2_x}\mathbf{\hat{a_2}_x}\otimes \mathbf{\hat{a_2}_x} + I_{2_y}\mathbf{\hat{a_2}_y}\otimes \mathbf{\hat{a_2}_y} + I_{2_z}\mathbf{\hat{a_2}_z}\otimes \mathbf{\hat{a_2}_z}



:math:`\text{Gyrostat}\ G_1`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow {F_G}^* = -m_G {a^G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_G1 = -M1 * N_a_Q
    Fstar_G1




.. math::

    \displaystyle -  M_{1} \dot{u}_{2}\mathbf{\hat{a_1}_x} -  M_{1} \left(u_{1} + \frac{u_{2} \operatorname{sin}\left(q_{4}\right)}{L}\right) u_{2}\mathbf{\hat{a_1}_y}



:math:`\rightarrow {T_G}^* \overset{\Delta}{=} -[\alpha_A \cdot I_G + {\omega_r}^A \times (I_G \cdot {\omega_r}^A)]`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Tstar_G1 = -(N_aa_A1.dot(I1) + me.cross(N_w_A1, I1.dot(N_w_A1)))
    Tstar_G1




.. math::

    \displaystyle -  I_{1_z} \left(\dot{u}_{1} + \frac{u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{\operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)\mathbf{\hat{a_1}_z}



:math:`\rightarrow ({F_r}^*)_{GR} = {V_r}^G \cdot {F_G}^* + {\omega_r}^A \cdot {T_G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_G1_R = V_1_Q.dot(Fstar_G1) + w_1_A1.dot(Tstar_G1).subs(soldict)
    Fstar_1_G1_R.subs({N_w_A1.dt(N).subs(soldict).dot(A1.z): alpha__A1})




.. math::

    \displaystyle - I_{1_z} \alpha^{A1}



.. code:: ipython3

    Fstar_2_G1_R = V_2_Q.dot(Fstar_G1) + w_2_A1.dot(Tstar_G1).subs(soldict)
    Fstar_2_G1_R.subs({N_w_A1.dt(N).subs(soldict).dot(A1.z): alpha__A1})




.. math::

    \displaystyle - \frac{I_{1_z} \alpha^{A1} \operatorname{sin}\left(q_{4}\right)}{L} - M_{1} \dot{u}_{2}



:math:`\rightarrow (F_r^*)_{GI} = -J\{\omega_r^A \cdot [\ddot{q_k} g_1 + \dot{q_k} (\omega_3^A g_2 - \omega_2^A g_3)] + C_{kr} (\dot{\omega}_1^A + \ddot{q_k}) \} \qquad (r=1,...,n-m)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`\text{Here}, \{\omega_1^A: \omega_2^A,\ \omega_2^A: \omega_3^A,\ \omega_3^A: \omega_1^A\}`

:math:`\rightarrow \dot{q_k} = \sum_{s = 1}^{n - m} C_{ks} u_s + D_k \quad \text{(Generalized Speeds)}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`\omega_i^A \overset{\Delta}{=} \omega^A \cdot \hat{g}_i \quad (i = 1, 2, 3)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # C_kr
    C51, C61 = sm.symbols('C_51, C_61')
    C_51 = soldict[q5.diff()].diff(u1)
    C_61 = soldict[q6.diff()].diff(u1)
    Fstar_1_G1_I = -J1 * (N_w_A1.dot(q5.diff().diff() * A1.y + q5.diff()*(N_w_A1.dot(A1.x)*A1.z - N_w_A1.dot(A1.z)*A1.x)) + C_51 * (N_w_A1.dot(A1.y).diff() + q5.diff().diff())) \
                   -J1 * (N_w_A1.dot(q6.diff().diff() * A1.y + q6.diff()*(N_w_A1.dot(A1.x)*A1.z - N_w_A1.dot(A1.z)*A1.x)) + C_61 * (N_w_A1.dot(A1.y).diff() + q6.diff().diff()))   # B1 \ B2
    
    Fstar_1_G1_I, C_51, C_61, Fstar_1_G1_I.subs({-C_51: -C51, -C_61: -C61}).simplify()




.. math::

    \displaystyle \left( \frac{J_{1} a \ddot{q}_{5}}{r_{1}} - \frac{J_{1} a \ddot{q}_{6}}{r_{1}}, \  - \frac{a}{r_{1}}, \  \frac{a}{r_{1}}, \  - J_{1} \left(C_{51} \ddot{q}_{5} + C_{61} \ddot{q}_{6}\right)\right)



.. code:: ipython3

    # C_kr 
    C52, C62 = sm.symbols('C_52, C_62')
    C_52 = soldict[q5.diff()].diff(u2)
    C_62 = soldict[q6.diff()].diff(u2)
    Fstar_2_G1_I = -J1 * (N_w_A1.dot(q5.diff().diff() * A1.y + q5.diff()*(N_w_A1.dot(A1.x)*A1.z - N_w_A1.dot(A1.z)*A1.x)) + C_52 * (N_w_A1.dot(A1.y).diff() + q5.diff().diff())) \
                   -J1 * (N_w_A1.dot(q6.diff().diff() * A1.y + q6.diff()*(N_w_A1.dot(A1.x)*A1.z - N_w_A1.dot(A1.z)*A1.x)) + C_62 * (N_w_A1.dot(A1.y).diff() + q6.diff().diff()))   # B1 \ B2
    
    Fstar_2_G1_I, C_52, C_62, Fstar_2_G1_I.subs({-C_52: -C52, -C_62: -C62}).simplify()




.. math::

    \displaystyle \left( \frac{J_{1} \left(- L + a \operatorname{sin}\left(q_{4}\right)\right) \ddot{q}_{5}}{L r_{1}} - \frac{J_{1} \left(L + a \operatorname{sin}\left(q_{4}\right)\right) \ddot{q}_{6}}{L r_{1}}, \  - \frac{- L + a \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  \frac{L + a \operatorname{sin}\left(q_{4}\right)}{L r_{1}}, \  - J_{1} \left(C_{52} \ddot{q}_{5} + C_{62} \ddot{q}_{6}\right)\right)



:math:`\rightarrow (F_r^*)_G = (F_r^*)_{GR} + (F_r^*)_{GI}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_G1 = Fstar_1_G1_R + Fstar_1_G1_I
    Fstar_1_G1.subs({N_w_A1.dt(N).subs(soldict).dot(A1.z): alpha__A1}).subs({-C_51: -C51, -C_61: -C61}).simplify()




.. math::

    \displaystyle - C_{51} J_{1} \ddot{q}_{5} - C_{61} J_{1} \ddot{q}_{6} - I_{1_z} \alpha^{A1}



.. code:: ipython3

    Fstar_2_G1 = Fstar_2_G1_R + Fstar_2_G1_I
    Fstar_2_G1.subs({N_w_A1.dt(N).subs(soldict).dot(A1.z): alpha__A1}).subs({-C_52: -C52, -C_62: -C62}).simplify()




.. math::

    \displaystyle - C_{52} J_{1} \ddot{q}_{5} - C_{62} J_{1} \ddot{q}_{6} - \frac{I_{1_z} \alpha^{A1} \operatorname{sin}\left(q_{4}\right)}{L} - M_{1} \dot{u}_{2}



:math:`\text{Gyrostat}\ G_2`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow {F_G}^* = -m_G {a^G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_G2 = -M2 * N_a_C
    Fstar_G2




.. math::

    \displaystyle -  M_{2} \left(- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2} - \frac{l u^{2}_{2} \operatorname{sin}^{2}\left(q_{4}\right)}{L^{2}}\right)\mathbf{\hat{a_2}_x} -  M_{2} \left(\frac{l u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{l \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L} + \frac{u^{2}_{2} \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{4}\right)}{L}\right)\mathbf{\hat{a_2}_y}



:math:`\rightarrow {T_G}^* \overset{\Delta}{=} -[\alpha_A \cdot I_G + {\omega_r}^A \times (I_G \cdot {\omega_r}^A)]`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Tstar_G2 = -(N_aa_A2.dot(I2) + me.cross(N_w_A2, I2.dot(N_w_A2)))
    Tstar_G2




.. math::

    \displaystyle -  I_{2_z} \left(\frac{u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{\operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)\mathbf{\hat{a_2}_z}



:math:`\rightarrow ({F_r}^*)_{GR} = {V_r}^G \cdot {F_G}^* + {\omega_r}^A \cdot {T_G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_G2_R = V_1_C.dot(Fstar_G2) + w_1_A2.dot(Tstar_G2).subs(soldict)
    Fstar_1_G2_R.subs({N_w_A2.dt(N).subs(soldict).dot(A2.z): alpha__A2})




.. math::

    \displaystyle 0



.. code:: ipython3

    Fstar_2_G2_R = V_2_C.dot(Fstar_G2) + w_2_A1.dot(Tstar_G2).subs(soldict)
    Fstar_2_G2_R.subs({N_w_A2.dt(N).subs(soldict).dot(A2.z): alpha__A2})




.. math::

    \displaystyle - \frac{I_{2_z} \alpha^{A2} \operatorname{sin}\left(q_{4}\right)}{L} - M_{2} \left(- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2} - \frac{l u^{2}_{2} \operatorname{sin}^{2}\left(q_{4}\right)}{L^{2}}\right) \operatorname{cos}\left(q_{4}\right) - \frac{M_{2} l \left(\frac{l u_{1} u_{2} \operatorname{cos}\left(q_{4}\right)}{L} + \frac{l \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L} + \frac{u^{2}_{2} \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{4}\right)}{L}\right) \operatorname{sin}\left(q_{4}\right)}{L}



:math:`\rightarrow (F_r^*)_{GI} = -J\{\omega_r^A \cdot [\ddot{q_k} g_1 + \dot{q_k} (\omega_3^A g_2 - \omega_2^A g_3)] + C_{kr} (\dot{\omega}_1^A + \ddot{q_k}) \} \qquad (r=1,...,n-m)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`\text{Here}, \{\omega_1^A: \omega_2^A,\ \omega_2^A: \omega_3^A,\ \omega_3^A: \omega_1^A\}`

:math:`\rightarrow \dot{q_k} = \sum_{s = 1}^{n - m} C_{ks} u_s + D_k \quad \text{(Generalized Speeds)}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`\omega_i^A \overset{\Delta}{=} \omega^A \cdot \hat{g}_i \quad (i = 1, 2, 3)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    # C_kr
    C71, C81 = sm.symbols('C_71, C_81')
    C_71 = soldict[q7.diff()].diff(u1)
    C_81 = soldict[q8.diff()].diff(u1)
    Fstar_1_G2_I = -J2 * (N_w_A2.dot(q7.diff().diff() * A2.y + q7.diff()*(N_w_A2.dot(A2.x)*A2.z - N_w_A2.dot(A2.z)*A2.x)) + C_71 * (N_w_A2.dot(A2.y).diff() + q7.diff().diff())) \
                   -J2 * (N_w_A2.dot(q8.diff().diff() * A2.y + q8.diff()*(N_w_A2.dot(A2.x)*A2.z - N_w_A2.dot(A2.z)*A2.x)) + C_81 * (N_w_A2.dot(A2.y).diff() + q8.diff().diff()))   # B1 \ B2
    
    Fstar_1_G2_I, C_71, C_81, # Fstar_1_G2_I.subs({-C_71: -C71, -C_81: -C81}).simplify()




.. math::

    \displaystyle \left( 0, \  0, \  0\right)



.. code:: ipython3

    # C_kr 
    C72, C82 = sm.symbols('C_72, C_82')
    C_72 = soldict[q7.diff()].diff(u2)
    C_82 = soldict[q8.diff()].diff(u2)
    Fstar_2_G2_I = -J2 * (N_w_A2.dot(q7.diff().diff() * A2.y + q7.diff()*(N_w_A2.dot(A2.x)*A2.z - N_w_A2.dot(A2.z)*A2.x)) + C_72 * (N_w_A2.dot(A2.y).diff() + q7.diff().diff())) \
                   -J2 * (N_w_A2.dot(q8.diff().diff() * A2.y + q8.diff()*(N_w_A2.dot(A2.x)*A2.z - N_w_A2.dot(A2.z)*A2.x)) + C_82 * (N_w_A2.dot(A2.y).diff() + q8.diff().diff()))   # B1 \ B2
    
    Fstar_2_G2_I, C_72, C_82, Fstar_2_G2_I.subs({-C_72: -C72, -C_82: -C82}).simplify()




.. math::

    \displaystyle \left( \frac{J_{2} \left(- L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) \ddot{q}_{7}}{L r_{2}} - \frac{J_{2} \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) \ddot{q}_{8}}{L r_{2}}, \  - \frac{- L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)}{L r_{2}}, \  \frac{L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)}{L r_{2}}, \  - J_{2} \left(C_{72} \ddot{q}_{7} + C_{82} \ddot{q}_{8}\right)\right)



:math:`\rightarrow (F_r^*)_G = (F_r^*)_{GR} + (F_r^*)_{GI}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_G2 = Fstar_1_G2_R + Fstar_1_G2_I
    # Fstar_1_G2.subs({N_w_A2.dt(N).subs(soldict).dot(A2.z): alpha__A2}) # .subs({-C_71: -C71, -C_81: -C81}).simplify()
    Fstar_1_G2 = 0

:math:`\text{Here}, \{a_1^C: a_2^C,\ a_2^C: a_3^C,\ a_3^C: a_1^C\}`

.. code:: ipython3

    Fstar_2_G2 = Fstar_2_G2_R + Fstar_2_G2_I
    Fstar_2_G2.subs({N_w_A2.dt(N).subs(soldict).dot(A2.z): alpha__A2}).subs({N_v_C.dt(N).subs(soldict).dot(A2.x): a_3__C}).subs({N_v_C.dt(N).subs(soldict).dot(A2.y): a_1__C}).subs({-C_72: -C72, -C_82: -C82}).simplify()




.. math::

    \displaystyle - C_{72} J_{2} \ddot{q}_{7} - C_{82} J_{2} \ddot{q}_{8} - \frac{I_{2_z} \alpha^{A2} \operatorname{sin}\left(q_{4}\right)}{L} - M_{2} a^{C}_{3} \operatorname{cos}\left(q_{4}\right) - \frac{M_{2} a^{C}_{1} l \operatorname{sin}\left(q_{4}\right)}{L}



:math:`\text{Variable-Mass Particle}\ P`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow {F_G}^* = -m_G {a^G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_P = -m * N_a_P
    Fstar_P




.. math::

    \displaystyle -  \left(- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2}\right) m\mathbf{\hat{a_2}_x} -  \frac{m u^{2}_{2} \operatorname{sin}\left(q_{4}\right) \operatorname{cos}\left(q_{4}\right)}{L}\mathbf{\hat{a_2}_y}



:math:`\rightarrow ({F_r}^*)_{GR} = {V_r}^G \cdot {F_G}^*`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_P_R = V_1_P.dot(Fstar_P)
    Fstar_1_P_R




.. math::

    \displaystyle 0



.. code:: ipython3

    Fstar_2_P_R = V_2_P.dot(Fstar_P) 
    Fstar_2_P_R




.. math::

    \displaystyle - \left(- u_{1} u_{2} \operatorname{sin}\left(q_{4}\right) + \operatorname{cos}\left(q_{4}\right) \dot{u}_{2}\right) m \operatorname{cos}\left(q_{4}\right)



:math:`\rightarrow (F_r^*)_G = (F_r^*)_{GR}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    Fstar_1_P = Fstar_1_P_R
    Fstar_1_P




.. math::

    \displaystyle 0



:math:`\text{Here}, \{a_1^P: a_2^P,\ a_2^P: a_3^P,\ a_3^P: a_1^P\}`

.. code:: ipython3

    Fstar_2_P = Fstar_2_P_R
    Fstar_2_P.subs({N_v_P.dt(N).subs(soldict).dot(A2.x): a_3__P}).subs({N_v_P.dt(N).subs(soldict).dot(A2.y): a_1__P}).simplify()




.. math::

    \displaystyle - a^{P}_{3} m \operatorname{cos}\left(q_{4}\right)



:math:`\text{Generalized Inertia Forces}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow F_r^* = (F_r^*)_{G_1} + (F_r^*)_{G_2} + (F_r^*)_{P} \quad (r = 1, 2)`

.. code:: ipython3

    Fstar_1 = Fstar_1_G1 + Fstar_1_G2 + Fstar_1_P
    Fstar_1.subs(soldict).simplify()




.. math::

    \displaystyle \frac{- I_{1_z} r_{1} \left(L \dot{u}_{1} + u_{1} u_{2} \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}\right) + J_{1} L a \left(\frac{- a \dot{u}_{1} + \dot{u}_{2} - \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} - \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}}{r_{1}} - \frac{a \dot{u}_{1} + \dot{u}_{2} + \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} + \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}}{r_{1}}\right)}{L r_{1}}



.. code:: ipython3

    Fstar_2 = Fstar_2_G1 + Fstar_2_G2 + Fstar_2_P
    Fstar_2.subs(soldict).simplify()




.. math::

    \displaystyle \frac{- J_{1} L r_{2} \left(\frac{\left(L - a \operatorname{sin}\left(q_{4}\right)\right) \left(- a \dot{u}_{1} + \dot{u}_{2} - \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} - \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)}{r_{1}} + \frac{\left(L + a \operatorname{sin}\left(q_{4}\right)\right) \left(a \dot{u}_{1} + \dot{u}_{2} + \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} + \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)}{r_{1}}\right) - J_{2} L r_{1} \left(\left(L \operatorname{cos}\left(q_{4}\right) - b \operatorname{sin}\left(q_{4}\right)\right) \left(\frac{\left(- \operatorname{sin}\left(q_{4}\right) \dot{q}_{4} - \frac{b \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L}\right) u_{2}}{r_{2}} + \frac{\left(\operatorname{cos}\left(q_{4}\right) - \frac{b \operatorname{sin}\left(q_{4}\right)}{L}\right) \dot{u}_{2}}{r_{2}}\right) + \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) \left(\frac{\left(- \operatorname{sin}\left(q_{4}\right) \dot{q}_{4} + \frac{b \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L}\right) u_{2}}{r_{2}} + \frac{\left(\operatorname{cos}\left(q_{4}\right) + \frac{b \operatorname{sin}\left(q_{4}\right)}{L}\right) \dot{u}_{2}}{r_{2}}\right)\right) + \frac{L^{2} r_{1} r_{2} \left(- 2 M_{1} \dot{u}_{2} + \left(u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) - \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - \dot{u}_{2}\right) m\right)}{2} + \frac{r_{1} r_{2} \left(- 2 I_{1_z} L \operatorname{sin}\left(q_{4}\right) \dot{u}_{1} - I_{1_z} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + I_{1_z} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - I_{1_z} \dot{u}_{2} - I_{2_z} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + I_{2_z} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - I_{2_z} \dot{u}_{2} + L^{2} M_{2} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) - L^{2} M_{2} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - L^{2} M_{2} \dot{u}_{2} - M_{2} l^{2} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + M_{2} l^{2} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - M_{2} l^{2} \dot{u}_{2}\right)}{2}}{L^{2} r_{1} r_{2}}



:math:`\mathbf{\text{Velocity of material ejected at}\ P\ \text{relative to}\ A_2 \rightarrow -C(t)g_3^{'};\ C(t) \rightarrow \text{positive}}`

.. code:: ipython3

    C = me.dynamicsymbols('C')
    C_t = -C*A2.x
    C_t




.. math::

    \displaystyle -  C\mathbf{\hat{a_2}_x}



:math:`\text{Generalized Thrust}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow F_r^{\prime} \triangleq \sum_{i=1}^{N} \mathbf{V}_{r}^{P i} \cdot \mathbf{C}^{P i} \dot{m}_{i} \quad (r=1, \ldots, k)`

.. code:: ipython3

    Fprime_1 = V_1_P.dot(C_t)*m.diff()
    Fprime_1




.. math::

    \displaystyle 0



.. code:: ipython3

    Fprime_2 = V_2_P.dot(C_t)*m.diff()
    Fprime_2




.. math::

    \displaystyle - C \operatorname{cos}\left(q_{4}\right) \dot{m}



:math:`\text{Extended Kane's Equations for Variable-Mass Systems}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\rightarrow F_r + F_r^* + F_r^{\prime} = 0 \quad (r = 1,..., k)`

Here, :math:`F_r = 0` :math:`\rightarrow` no forces contributing to
generalized active forces

.. code:: ipython3

    kane_1 = Fstar_1.simplify() + Fprime_1.simplify()
    kane_1.subs(soldict).simplify()




.. math::

    \displaystyle \frac{- I_{1_z} r_{1} \left(L \dot{u}_{1} + u_{1} u_{2} \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}\right) + J_{1} L a \left(\frac{- a \dot{u}_{1} + \dot{u}_{2} - \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} - \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}}{r_{1}} - \frac{a \dot{u}_{1} + \dot{u}_{2} + \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} + \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}}{r_{1}}\right)}{L r_{1}}



.. code:: ipython3

    kane_2 = Fstar_2 + Fprime_2
    kane_2.subs(soldict).simplify()




.. math::

    \displaystyle \frac{- J_{1} L r_{2} \left(\frac{\left(L - a \operatorname{sin}\left(q_{4}\right)\right) \left(- a \dot{u}_{1} + \dot{u}_{2} - \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} - \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)}{r_{1}} + \frac{\left(L + a \operatorname{sin}\left(q_{4}\right)\right) \left(a \dot{u}_{1} + \dot{u}_{2} + \frac{a u_{2} \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L} + \frac{a \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}}{L}\right)}{r_{1}}\right) - J_{2} L r_{1} \left(\left(L \operatorname{cos}\left(q_{4}\right) - b \operatorname{sin}\left(q_{4}\right)\right) \left(\frac{\left(- \operatorname{sin}\left(q_{4}\right) \dot{q}_{4} - \frac{b \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L}\right) u_{2}}{r_{2}} + \frac{\left(\operatorname{cos}\left(q_{4}\right) - \frac{b \operatorname{sin}\left(q_{4}\right)}{L}\right) \dot{u}_{2}}{r_{2}}\right) + \left(L \operatorname{cos}\left(q_{4}\right) + b \operatorname{sin}\left(q_{4}\right)\right) \left(\frac{\left(- \operatorname{sin}\left(q_{4}\right) \dot{q}_{4} + \frac{b \operatorname{cos}\left(q_{4}\right) \dot{q}_{4}}{L}\right) u_{2}}{r_{2}} + \frac{\left(\operatorname{cos}\left(q_{4}\right) + \frac{b \operatorname{sin}\left(q_{4}\right)}{L}\right) \dot{u}_{2}}{r_{2}}\right)\right) + \frac{L^{2} r_{1} r_{2} \left(- 2 M_{1} \dot{u}_{2} + \left(u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) - \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - \dot{u}_{2}\right) m - 2 C \operatorname{cos}\left(q_{4}\right) \dot{m}\right)}{2} + \frac{r_{1} r_{2} \left(- 2 I_{1_z} L \operatorname{sin}\left(q_{4}\right) \dot{u}_{1} - I_{1_z} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + I_{1_z} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - I_{1_z} \dot{u}_{2} - I_{2_z} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + I_{2_z} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - I_{2_z} \dot{u}_{2} + L^{2} M_{2} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) - L^{2} M_{2} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - L^{2} M_{2} \dot{u}_{2} - M_{2} l^{2} u_{1} u_{2} \operatorname{sin}\left(2 q_{4}\right) + M_{2} l^{2} \operatorname{cos}\left(2 q_{4}\right) \dot{u}_{2} - M_{2} l^{2} \dot{u}_{2}\right)}{2}}{L^{2} r_{1} r_{2}}



.. code:: ipython3

    kane_1_eq = sm.Eq(kane_1.simplify().subs(soldict).simplify().subs(u).simplify(), 0)
    kane_1_eq




.. math::

    \displaystyle - \frac{\left(I_{1_z} r_{1}^{2} + 2 J_{1} a^{2}\right) \left(L \dot{u}_{1} + u_{1} u_{2} \operatorname{cos}\left(q_{4}\right) + \operatorname{sin}\left(q_{4}\right) \dot{u}_{2}\right)}{L r_{1}^{2}} = 0



.. code:: ipython3

    kane_2_eq = sm.Eq(kane_2.simplify().subs(soldict).simplify().subs(u).simplify(), 0)
    # kane_2_eq
