==================================================
Exercises from Chapter 2 in Kane and Levinson 1985
==================================================

Exercise 2.7
============

.. jupyter-kernel:: python3
   :id: ex2-7

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex2-7` or Jupyter notebook:
   :jupyter-download:notebook:`ex2-7`.

.. jupyter-execute::

   from sympy.physics.mechanics import ReferenceFrame, dynamicsymbols, mprint
   from sympy import solve, pi, Eq

   q1, q2, q3, q4, q5 = dynamicsymbols('q1 q2 q3 q4 q5')
   q1d, q2d, q3d, q4d, q5d = dynamicsymbols('q1 q2 q3 q4 q5', level=1)

   ux, uy, uz = dynamicsymbols('ux uy uz')
   u1, u2, u3 = dynamicsymbols('u1 u2 u3')

   A = ReferenceFrame('A')
   B_prime = A.orientnew('B_prime', 'Axis', [q1, A.z])
   B = B_prime.orientnew('B', 'Axis', [pi/2 - q2, B_prime.x])
   C = B.orientnew('C', 'Axis', [q3, B.z])

   # Angular velocity based on coordinate time derivatives
   w_C_in_A_qd = C.ang_vel_in(A)

   # First definition of Angular velocity
   w_C_in_A_uxuyuz = ux * A.x + uy * A.y + uz * A.z
   print("Using w_C_A as")
   print(w_C_in_A_uxuyuz)

.. jupyter-execute::

   kinematic_eqs = [(w_C_in_A_qd - w_C_in_A_uxuyuz) & uv for uv in A]
   print("The kinematic equations are:")
   soln = solve(kinematic_eqs, [q1d, q2d, q3d])
   for qd in [q1d, q2d, q3d]:
       mprint(Eq(qd, soln[qd]))

.. jupyter-execute::

   # Second definition of Angular velocity
   w_C_in_A_u1u2u3 = u1 * B.x + u2 * B.y + u3 * B.z
   print("Using w_C_A as")
   print(w_C_in_A_u1u2u3)

.. jupyter-execute::

   kinematic_eqs = [(w_C_in_A_qd - w_C_in_A_u1u2u3) & uv for uv in A]
   print("The kinematic equations are:")
   soln = solve(kinematic_eqs, [q1d, q2d, q3d])
   for qd in [q1d, q2d, q3d]:
       mprint(Eq(qd, soln[qd]))

Exercise 3.10
=============

.. jupyter-kernel:: python3
   :id: ex3-10

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex3-10` or Jupyter notebook:
   :jupyter-download:notebook:`ex3-10`.

.. jupyter-execute::

   from sympy import cancel, collect, expand_trig, solve, symbols, trigsimp
   from sympy import sin, cos
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import dot, dynamicsymbols, msprint

   q1, q2, q3, q4, q5, q6, q7 = q = dynamicsymbols('q1:8')
   u1, u2, u3, u4, u5, u6, u7 = u = dynamicsymbols('q1:8', level=1)

   r, theta, b = symbols('r θ b', real=True, positive=True)

   # define reference frames
   R = ReferenceFrame('R') # fixed race rf, let R.z point upwards
   A = R.orientnew('A', 'axis', [q7, R.z]) # rf that rotates with S* about R.z
   # B.x, B.z are parallel with face of cone, B.y is perpendicular
   B = A.orientnew('B', 'axis', [-theta, A.x])
   S = ReferenceFrame('S')
   S.set_ang_vel(A, u1*A.x + u2*A.y + u3*A.z)
   C = ReferenceFrame('C')
   C.set_ang_vel(A, u4*B.x + u5*B.y + u6*B.z)

   # define points
   pO = Point('O')
   pS_star = pO.locatenew('S*', b*A.y)
   pS_hat = pS_star.locatenew('S^', -r*B.y) # S^ touches the cone
   pS1 = pS_star.locatenew('S1', -r*A.z) # S1 touches horizontal wall of the race
   pS2 = pS_star.locatenew('S2', r*A.y) # S2 touches vertical wall of the race

   pO.set_vel(R, 0)
   pS_star.v2pt_theory(pO, R, A)
   pS1.v2pt_theory(pS_star, R, S)
   pS2.v2pt_theory(pS_star, R, S)

   # Since S is rolling against R, v_S1_R = 0, v_S2_R = 0.
   vc = [dot(p.vel(R), basis) for p in [pS1, pS2] for basis in R]

   pO.set_vel(C, 0)
   pS_star.v2pt_theory(pO, C, A)
   pS_hat.v2pt_theory(pS_star, C, S)

   # Since S is rolling against C, v_S^_C = 0.
   # Cone has only angular velocity in R.z direction.
   vc += [dot(pS_hat.vel(C), basis) for basis in A]
   vc += [dot(C.ang_vel_in(R), basis) for basis in [R.x, R.y]]
   vc_map = solve(vc, u)

   # Pure rolling between S and C, dot(ω_C_S, B.y) = 0.
   b_val = solve([dot(C.ang_vel_in(S), B.y).subs(vc_map).simplify()], b)[0][0]
   print('b = {0}'.format(msprint(collect(cancel(expand_trig(b_val)), r))))

.. jupyter-execute::

   b_expected = r*(1 + sin(theta))/(cos(theta) - sin(theta))
   assert trigsimp(b_val - b_expected) == 0

Exercise 3.15
=============

.. jupyter-kernel:: python3
   :id: ex3-15

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex3-15` or Jupyter notebook:
   :jupyter-download:notebook:`ex3-15`.

A robot arm, composed of Rigid Bodies 'A', 'B', 'C', operates in Reference
Frame E. 'A*', 'B*', 'C*' are Points marking the centers of mass for the Rigid
Bodies 'A', 'B', 'C'.

Rigid Body 'D' also is lies in Reference Frame 'E'. The center of mass of 'D'
is marked as Point 'D*'. 'D' is fixed relative to 'C'.

Each Reference Frame has a set of mutually perpendicular vectors x, y, z.  'A'
is rotated by 'q0' relative to 'E' about an axis parallel to A.x. 'B' is
rotated by 'q1' relative to 'A' about an axis parallel to A.y. Point 'P' is
fixed in both 'A' and 'B'. A.x is parallel to E.x. A.y is parallel to B.y.

Point 'O' is a point fixed in both 'E' and 'A'. 'A*' is separated from 'O' by
'LA' * A.z. 'P' is separated from 'O' by 'LP' * A.z. 'B*' is separated from 'P'
by 'LB' * B.z. 'C*' is separated from 'B*' by 'q2' * B.z. 'D*' is separated
from 'C*' by ``p1*B.x + p2*B.y + p3*B.z``.

We define:'q0d' = 'u1', 'q1d' = 'u2', 'q2d' = 'u3'.  'LA', 'LB', 'LP', 'p1',
'p2', 'p3' are constants.

.. jupyter-execute::

   from sympy.physics.mechanics import dynamicsymbols, msprint
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy import solve, symbols

   # Define generalized coordinates, speeds, and constants:
   q0, q1, q2 = dynamicsymbols('q0 q1 q2')
   q0d, q1d, q2d = dynamicsymbols('q0 q1 q2', level=1)
   u1, u2, u3 = dynamicsymbols('u1 u2 u3')
   LA, LB, LP = symbols('LA LB LP')
   p1, p2, p3 = symbols('p1 p2 p3')

   E = ReferenceFrame('E')
   # A.x of Rigid Body A is fixed in Reference Frame E and is rotated by q0.
   A = E.orientnew('A', 'Axis', [q0, E.x])
   # B.y of Rigid Body B is fixed in Reference Frame A and is rotated by q1.
   B = A.orientnew('B', 'Axis', [q1, A.y])
   # Reference Frame C has no rotation relative to Reference Frame B.
   C = B.orientnew('C', 'Axis', [0, B.x])
   # Reference Frame D has no rotation relative to Reference Frame C.
   D = C.orientnew('D', 'Axis', [0, C.x])

   pO = Point('O')
   # The vector from Point O to Point A*, the center of mass of A, is LA * A.z.
   pAs = pO.locatenew('A*', LA * A.z)
   # The vector from Point O to Point P, which lies on the axis where
   # B rotates about A, is LP * A.z.
   pP = pO.locatenew('P', LP * A.z)
   # The vector from Point P to Point B*, the center of mass of B, is LB * B.z.
   pBs = pP.locatenew('B*', LB * B.z)
   # The vector from Point B* to Point C*, the center of mass of C, is q2 * B.z.
   pCs = pBs.locatenew('C*', q2 * B.z)
   # The vector from Point C* to Point D*, the center of mass of D,
   # is p1 * B.x + p2 * B.y + p3 * B.z.
   pDs = pCs.locatenew('D*', p1 * B.x + p2 * B.y + p3 * B.z)

   # Define generalized speeds as:
   # u1 = q0d
   # u2 = q1d
   # u3 = q2d
   A.set_ang_vel(E, u1 * A.x) # A.x = E.x
   B.set_ang_vel(A, u2 * B.y) # B.y = A.y
   pCs.set_vel(B, u3 * B.z)

   pO.set_vel(E, 0) # Point O is fixed in Reference Frame E
   pAs.v2pt_theory(pO, E, A) # Point A* is fixed in Reference Frame A
   pP.v2pt_theory(pO, E, A) # Point P is fixed in Reference Frame A
   pBs.v2pt_theory(pP, E, B) # Point B* is fixed in Reference Frame B
   pCs.v1pt_theory(pBs, E, B) # Point C* is moving in Reference Frame B
   pDs.set_vel(B, pCs.vel(B)) # Point D* is fixed relative to Point C* in B
   pDs.v1pt_theory(pBs, E, B) # Point D* is moving in Reference Frame B

   # Write generalized speeds as kinematic equations:
   kinematic_eqs = []
   kinematic_eqs.append(u1 - q0d)
   kinematic_eqs.append(u2 - q1d)
   kinematic_eqs.append(u3 - q2d)
   soln = solve(kinematic_eqs, [q0d, q1d, q2d])
   print("kinematic equations:")
   for qd in [q0d, q1d, q2d]:
      print("{0} = {1}".format(msprint(qd), msprint(soln[qd])))

.. jupyter-execute::

   ang_vels = ["\nangular velocities:"]
   ang_accs = ["\nangular accelerations:"]
   for rf in [A, B, C, D]:
      ang_v = getattr(rf, 'ang_vel_in')(E)
      ang_a = getattr(rf, 'ang_acc_in')(E)
      express_rf = B
      if rf == A:
         express_rf = A
      ang_vels.append("ang vel {0} wrt {1} = {2}".format(
               rf, E, ang_v.express(express_rf)))
      ang_accs.append("ang acc {0} wrt {1} = {2}".format(
               rf, E, ang_a.express(express_rf)))

   vels = ["\nvelocities:"]
   accs = ["\naccelerations:"]
   for point in [pAs, pBs, pCs, pDs]:
      v = getattr(point, 'vel')(E)
      a = getattr(point, 'acc')(E)
      express_rf = B
      if point == pAs:
         express_rf = A
      vels.append("vel {0} wrt {1} = {2}".format(
               point, E, v.express(express_rf)))
      accs.append("acc {0} wrt {1} = {2}".format(
               point, E, a.express(express_rf)))

   for results in ang_vels + ang_accs + vels + accs:
      print(results)

Exercise 4.1
============

.. jupyter-kernel:: python3
   :id: ex4-1

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex4-1` or Jupyter notebook:
   :jupyter-download:notebook:`ex4-1`.

.. jupyter-execute::

   from sympy.physics.mechanics import dot, msprint
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy import symbols, sin, cos
   from sympy.simplify.simplify import trigsimp

   theta = symbols('theta:3')
   x = symbols('x:3')
   q = symbols('q')

   A = ReferenceFrame('A')
   B = A.orientnew('B', 'SPACE', theta, 'xyz')

   O = Point('O')
   P = O.locatenew('P', x[0] * A.x + x[1] * A.y + x[2] * A.z)
   p = P.pos_from(O)

   # From problem, point P is on L (span(B.x)) when:
   constraint_eqs = {x[0] : q*cos(theta[1])*cos(theta[2]),
                     x[1] : q*cos(theta[1])*sin(theta[2]),
                     x[2] : -q*sin(theta[1])}

   # If point P is on line L then r^{P/O} will have no components in the B.y or
   # B.z directions since point O is also on line L and B.x is parallel to L.
   assert(trigsimp(dot(P.pos_from(O), B.y).subs(constraint_eqs)) == 0)
   assert(trigsimp(dot(P.pos_from(O), B.z).subs(constraint_eqs)) == 0)

Exercise 4.18
=============

.. jupyter-kernel:: python3
   :id: ex4-18

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex4-18` or Jupyter notebook:
   :jupyter-download:notebook:`ex4-18`.

.. jupyter-execute::

   from sympy.physics.mechanics import dynamicsymbols, msprint
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy import solve, symbols, pi

   # Define generalized coordinates, speeds, and constants:
   qi = dynamicsymbols('q0 q1 q2 q3 q4 q5')
   qid = dynamicsymbols('q0 q1 q2 q3 q4 q5', level=1)
   ui = dynamicsymbols('u0 u1 u2 u3 u4 u5')
   R = symbols('R')

   A = ReferenceFrame('A')
   A_1 = A.orientnew("A'", 'Axis', [qi[1], A.z])
   B = A_1.orientnew('B', 'Axis', [pi/2 - qi[2], A_1.x])
   C = B.orientnew('C', 'Axis', [qi[3], B.z])

   pO = Point('O')
   pCs = pO.locatenew('C*', qi[4] * A.x + qi[5] * A.y + R * B.y)

   pO.set_vel(A, 0) # Point O is fixed in Reference Frame A
   pCs.v2pt_theory(pO, A, B) # Point C* is fixed in Reference Frame B

   # Set ui = qid
   kinematic_eqs = []
   for u, qd in zip(ui, qid):
      kinematic_eqs.append(u - qd)
   soln = solve(kinematic_eqs, qid)
   print("kinematic equations:")
   for qd in qid:
      print("{0} = {1}".format(msprint(qd), msprint(soln[qd])))

   print("\nposition of C* from O = {0}".format(msprint(pCs.pos_from(pO))))
   print("\nvelocity of C* wrt A = {0}".format(msprint(pCs.vel(A).express(B))))
