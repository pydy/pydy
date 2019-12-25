==================================================
Exercises from Chapter 3 in Kane and Levinson 1985
==================================================

Exercise 5.1
============

.. jupyter-kernel:: python3
   :id: ex5-1

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex5-1` or Jupyter notebook:
   :jupyter-download:notebook:`ex5-1`.

.. jupyter-execute::

   from sympy import symbols, integrate
   from sympy import pi, acos
   from sympy import Matrix, S, simplify
   from sympy.physics.mechanics import dot, ReferenceFrame
   import scipy.integrate


   def eval_v(v, N):
      return sum(dot(v, n).evalf() * n for n in N)


   def integrate_v(integrand, rf, bounds):
      """Return the integral for a Vector integrand."""
      # integration problems if option meijerg=True is not used
      return sum(simplify(integrate(dot(integrand, n),
                                    bounds,
                                    meijerg=True)) * n for n in rf)

   m, R = symbols('m R', real=True, nonnegative=True)
   theta, r, x, y, z = symbols('theta r x y z', real=True)

Regarding Fig. P5.1 as showing two views of a body B formed by matter
distributed uniformly (a) over a surface having no planar portions and (b)
throughout a solid, determine (by integration) the coordinates x*, y*, z* of
the mass center of B.

.. jupyter-execute::

   A = ReferenceFrame('A')

   # part a
   print("a) mass center of surface area")
   SA = (2 * pi * R) * 2
   rho_a = m / SA

   B = A.orientnew('B', 'Axis', [theta, A.x])
   p = x * A.x  + r * B.y
   y_ = dot(p, A.y)
   z_ = dot(p, A.z)
   J = Matrix([x, y_, z_]).jacobian([x, r, theta])
   dJ = simplify(J.det())
   print("dJ = {0}".format(dJ))

   # calculate center of mass for the cut cylinder
   # ranges from x = [1, 3], y = [-1, 1], z = [-1, 1] in Fig. P5.1
   mass_cc_a = rho_a * 2*pi*R
   cm_cc_a = (integrate_v(integrate_v((rho_a * p * dJ).subs(r, R),
                                      A, (theta, acos(1-x),
                                      2*pi - acos(1-x))),
                                      A, (x, 0, 2)) / mass_cc_a + A.x)
   print("cm = {0}".format(cm_cc_a))

.. jupyter-execute::

   mass_cyl_a = rho_a * 2*pi*R
   cm_cyl_a = A.x/S(2)

   cm_a = (mass_cyl_a*cm_cyl_a + mass_cc_a*cm_cc_a) / m
   print("center of mass = {0}".format(cm_a.subs(R, 1)))
   # part b
   print("b) mass center of volume")
   V = (pi*R**2) * 2
   rho_b = m / V
   mass_cc_b = rho_b * pi*R**2

   # calculate center of mass for the cut cylinder
   # ranges from x = [1, 3], y = [-1, 1], z = [-1, 1] in Fig. P5.1
   # compute the value using scipy due to issues with sympy
   R_ = 1
   def pos(z, y, x, coord):
      if coord == 0:
         return x
      elif coord == 1:
         return y
      elif coord == 2:
         return z
      else:
         raise ValueError
   # y bounds
   def ybu(x):
      return R_*(1 - x)
   def ybl(x):
      return -R_
   # z bounds
   def zbu(x, y):
      return (R_**2 - y**2)**0.5
   def zbl(x, y):
      return -1 * zbu(x, y)

   cm_cc_b = 0
   for i, b in enumerate(A):
      p_i = lambda z, y, x: pos(z, y, x, i)
      cm_cc_b += scipy.integrate.tplquad(p_i, 0, 2, ybl, ybu, zbl, zbu)[0] * b
   cm_cc_b *= (rho_b / mass_cc_b).subs(R, R_)
   cm_cc_b += A.x

   #integrand = rho_b * (x*A.x + y*A.y + z*A.z)
   #cm_cc_b = (integrate_v(integrate_v(integrate_v(integrand,
   #                                               A, (z,
   #                                                   -sqrt(R**2 - y**2),
   #                                                   sqrt(R**2 - y**2))),
   #                                   A, (y, -R, R*(1 - x))),
   #                       A, (x, 0, 2)) / mass_cc_b +
   #           A.x)
   print("cm = {0}".format(cm_cc_b))

   mass_cyl_b = rho_b * pi*R**2
   cm_cyl_b = A.x/S(2)

   cm_b = (mass_cyl_b*cm_cyl_b + mass_cc_b*cm_cc_b) / m
   print("center of mass = {0}".format(eval_v(cm_b.subs(R, 1), A)))

Exercise 5.4
============

.. jupyter-kernel:: python3
   :id: ex5-4

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex5-4` or Jupyter notebook:
   :jupyter-download:notebook:`ex5-4`.

.. jupyter-execute::

   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import dot
   from sympy import symbols, integrate
   from sympy import sin, cos, pi

   a, b, R = symbols('a b R', real=True, nonnegative=True)
   theta, x, y, z = symbols('theta x y z', real=True)
   ab_val = {a: 0.3, b: 0.3}

   # R is the radius of the circle, theta is the half angle.
   centroid_sector = 2*R*sin(theta) / (3 * theta)

   # common R, theta values
   theta_pi_4 = {theta: pi/4, R: a}
   R_theta_val = {theta: pi/4 * (1 - z/a), R: a}

   N = ReferenceFrame('N')
   def eval_vec(v):
       vs = v.subs(ab_val)
       return sum(dot(vs, n).evalf() * n for n in N)

   # For each part A, B, C, D, define an origin for that part such that the
   # centers of mass of each part of the component have positive N.x, N.y,
   # and N.z values.
   ## FIND CENTER OF MASS OF A
   vol_A_1 = pi * a**2 * b / 4
   vol_A_2 = pi * a**2 * a / 4 / 2
   vol_A = vol_A_1 + vol_A_2
   pA_O = Point('A_O')
   pAs_1 = pA_O.locatenew(
         'A*_1', (b/2 * N.z +
                  centroid_sector.subs(theta_pi_4) * sin(pi/4) * (N.x + N.y)))
   pAs_2 = pA_O.locatenew(
         'A*_2', (b * N.z +
                  N.z * integrate((theta*R**2*(z)).subs(R_theta_val),
                                    (z, 0, a)) / vol_A_2 +
                  N.x * integrate((theta*R**2 * cos(theta) *
                                    centroid_sector).subs(R_theta_val),
                                    (z, 0, a)) / vol_A_2 +
                  N.y * integrate((2*R**3/3 * 4*a/pi *
                                    sin(theta)**2).subs(R, a),
                                    (theta, 0, pi/4)) / vol_A_2))
   pAs = pA_O.locatenew('A*', ((pAs_1.pos_from(pA_O) * vol_A_1 +
                              pAs_2.pos_from(pA_O) * vol_A_2) /
                              vol_A))
   print('A* = {0}'.format(pAs.pos_from(pA_O)))
   print('A* = {0}'.format(eval_vec(pAs.pos_from(pA_O))))

   ## FIND CENTER OF MASS OF B
   vol_B_1 = pi*a**2/2
   vol_B_2 = a**2 / 2
   vol_B  = vol_B_1 + vol_B_2
   pB_O = Point('B_O')
   pBs_1 = pB_O.locatenew(
         'B*_1', (a*(N.x + N.z) + a/2*N.y +
                  (-N.x + N.z) * (R*sin(theta)/theta *
                                    sin(pi/4)).subs(theta_pi_4)))
   pBs_2 = pB_O.locatenew('B*_2', (a*N.y + a*N.z -
                                 (a/3 * N.y + a/3 * N.z)))
   pBs = pB_O.locatenew('B*', ((pBs_1.pos_from(pB_O) * vol_B_1 +
                              pBs_2.pos_from(pB_O) * vol_B_2) /
                              vol_B))
   print('\nB* = {0}'.format(pBs.pos_from(pB_O)))
   print('B* = {0}'.format(eval_vec(pBs.pos_from(pB_O))))

   ## FIND CENTER OF MASS OF C
   vol_C_1 = 2 * a**2 * b
   vol_C_2 = a**3 / 2
   vol_C_3 = a**3
   vol_C_4 = -pi*a**3/4
   vol_C = vol_C_1 + vol_C_2 + vol_C_3 + vol_C_4
   pC_O = Point('C_O')
   pCs_1 = pC_O.locatenew('C*_1', (a*N.x + a/2*N.y + b/2*N.z))
   pCs_2 = pC_O.locatenew('C*_2', (a*N.x + b*N.z +
                                 (a/3*N.x + a/2*N.y + a/3*N.z)))
   pCs_3 = pC_O.locatenew('C*_3', (b*N.z + a/2*(N.x + N.y + N.z)))
   pCs_4 = pC_O.locatenew(
         'C*_4', ((a + b)*N.z + a/2*N.y +
                  (N.x - N.z)*(centroid_sector.subs(
                           theta_pi_4)*sin(pi/4))))
   pCs = pC_O.locatenew('C*', ((pCs_1.pos_from(pC_O)*vol_C_1 +
                              pCs_2.pos_from(pC_O)*vol_C_2 +
                              pCs_3.pos_from(pC_O)*vol_C_3 +
                              pCs_4.pos_from(pC_O)*vol_C_4) /
                              vol_C))
   print('\nC* = {0}'.format(pCs.pos_from(pC_O)))
   print('C* = {0}'.format(eval_vec(pCs.pos_from(pC_O))))

   ## FIND CENTER OF MASS OF D
   vol_D = pi*a**3/4
   pD_O = Point('D_O')
   pDs = pD_O.locatenew('D*', (a*N.z + a/2*N.y +
                              (N.x - N.z)*(centroid_sector.subs(
                                       theta_pi_4) * sin(pi/4))))
   print('\nD* = {0}'.format(pDs.pos_from(pD_O)))
   print('D* = {0}'.format(eval_vec(pDs.pos_from(pD_O))))

   ## FIND CENTER OF MASS OF ASSEMBLY
   pO = Point('O')
   pA_O.set_pos(pO, 2*a*N.x - (a+b)*N.z)
   pB_O.set_pos(pO, 2*a*N.x - a*N.z)
   pC_O.set_pos(pO, -(a+b)*N.z)
   pD_O.set_pos(pO, -a*N.z)

   density_A = 7800
   density_B = 17.00
   density_C = 2700
   density_D = 8400
   mass_A = vol_A * density_A
   mass_B = vol_B * density_B
   mass_C = vol_C * density_C
   mass_D = vol_D * density_D

   pms = pO.locatenew('m*', ((pAs.pos_from(pO)*mass_A + pBs.pos_from(pO)*mass_B +
                              pCs.pos_from(pO)*mass_C + pDs.pos_from(pO)*mass_D) /
                           (mass_A + mass_B + mass_C + mass_D)))
   print('\nm* = {0}'.format(eval_vec(pms.pos_from(pO))))

Exercise 5.8
============

.. jupyter-kernel:: python3
   :id: ex5-8

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex5-8` or Jupyter notebook:
   :jupyter-download:notebook:`ex5-8`.

.. jupyter-execute::

   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import cross, dot
   from sympy import integrate, simplify, symbols, integrate
   from sympy import sin, cos, pi


   def calc_inertia_vec(rho, p, n_a, N, iv):
       integrand = rho * cross(p, cross(n_a, p))
       return sum(simplify(integrate(dot(integrand, n), iv)) * n
                  for n in N)


   a, b, L, l, m, h = symbols('a b L l m h', real=True, nonnegative=True)
   theta = symbols('theta', real=True)
   h_theta_val = {h:b*l/L, theta:2*pi*l/L}

   density = m/L
   N = ReferenceFrame('N')
   pO = Point('O')
   pP = pO.locatenew('P', h*N.x + a*cos(theta)*N.y + a*sin(theta)*N.z)

   I_1 = calc_inertia_vec(density, pP.pos_from(pO).subs(h_theta_val),
                          N.x, N, (l, 0, L))
   I_2 = calc_inertia_vec(density, pP.pos_from(pO).subs(h_theta_val),
                          N.y, N, (l, 0, L))
   print('I_1 = {0}'.format(I_1))
   print('I_2 = {0}'.format(I_2))

Exercise 5.12
=============

.. jupyter-kernel:: python3
   :id: ex5-12

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex5-12` or Jupyter notebook:
   :jupyter-download:notebook:`ex5-12`.

.. jupyter-execute::

   import numpy as np

   # n_a = 3/5*n_1 - 4/5*n_3. Substituting n_i for e_i results in
   # n_a = 4/5*e_1 + 3/5*n_2.
   a = np.matrix([4/5, 3/5, 0])
   Iij = np.matrix([[169, 144, -96],
                    [144, 260, 72],
                    [-96, 72, 325]])

   print("Moment of inertia of B with respect to a line that is parallel to")
   print("line PQ and passes through point O.")
   print("{0} kg m**2".format((a * Iij * a.T).item(0)))

Exercise 6.6
============

.. jupyter-kernel:: python3
   :id: ex6-6

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex6-6` or Jupyter notebook:
   :jupyter-download:notebook:`ex6-6`.

.. jupyter-execute::

   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import inertia, inertia_of_point_mass
   from sympy.physics.mechanics import dot
   from sympy import symbols
   from sympy import S

   m = symbols('m')
   m_val = 12

   N = ReferenceFrame('N')
   pO = Point('O')
   pBs = pO.locatenew('B*', -3*N.x + 2*N.y - 4*N.z)

   I_B_O = inertia(N, 260*m/m_val, 325*m/m_val, 169*m/m_val,
                   72*m/m_val, 96*m/m_val, -144*m/m_val)
   print("I_B_rel_O = {0}".format(I_B_O))

   I_Bs_O = inertia_of_point_mass(m, pBs.pos_from(pO), N)
   print("\nI_B*_rel_O = {0}".format(I_Bs_O))

   I_B_Bs = I_B_O - I_Bs_O
   print("\nI_B_rel_B* = {0}".format(I_B_Bs))

   pQ = pO.locatenew('Q', -4*N.z)
   I_Bs_Q = inertia_of_point_mass(m, pBs.pos_from(pQ), N)
   print("\nI_B*_rel_Q = {0}".format(I_Bs_Q))

   I_B_Q = I_B_Bs + I_Bs_Q
   print("\nI_B_rel_Q = {0}".format(I_B_Q))

   # n_a is a vector parallel to line PQ
   n_a = S(3)/5 * N.x - S(4)/5 * N.z
   I_a_a_B_Q = dot(dot(n_a, I_B_Q), n_a)
   print("\nn_a = {0}".format(n_a))
   print("\nI_a_a_B_Q = {0} = {1}".format(I_a_a_B_Q, I_a_a_B_Q.subs(m, m_val)))

Exercise 6.7
============

.. jupyter-kernel:: python3
   :id: ex6-7

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex6-7` or Jupyter notebook:
   :jupyter-download:notebook:`ex6-7`.

.. jupyter-execute::

   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import inertia, inertia_of_point_mass
   from sympy.physics.mechanics import dot
   from sympy import solve, sqrt, symbols, diff
   from sympy import sin, cos, pi, Matrix
   from sympy import simplify

   b, m, k_a = symbols('b m k_a', real=True, nonnegative=True)
   theta = symbols('theta', real=True)

   N = ReferenceFrame('N')
   N2 = N.orientnew('N2', 'Axis', [theta, N.z])
   pO = Point('O')
   pP1s = pO.locatenew('P1*', b/2 * (N.x + N.y))
   pP2s = pO.locatenew('P2*', b/2 * (2*N.x + N.y + N.z))
   pP3s = pO.locatenew('P3*', b/2 * ((2 + sin(theta))*N.x +
                                    (2 - cos(theta))*N.y +
                                    N.z))

   I_1s_O = inertia_of_point_mass(m, pP1s.pos_from(pO), N)
   I_2s_O = inertia_of_point_mass(m, pP2s.pos_from(pO), N)
   I_3s_O = inertia_of_point_mass(m, pP3s.pos_from(pO), N)
   print("\nI_1s_rel_O = {0}".format(I_1s_O))
   print("\nI_2s_rel_O = {0}".format(I_2s_O))
   print("\nI_3s_rel_O = {0}".format(I_3s_O))


   I_1_1s = inertia(N, m*b**2/12, m*b**2/12, 2*m*b**2/12)
   I_2_2s = inertia(N, 2*m*b**2/12, m*b**2/12, m*b**2/12)

   I_3_3s_N2 = Matrix([[2, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])
   I_3_3s_N = m*b**2/12 * simplify(N.dcm(N2) * I_3_3s_N2 * N2.dcm(N))
   I_3_3s = inertia(N, I_3_3s_N[0, 0], I_3_3s_N[1, 1], I_3_3s_N[2, 2],
                    I_3_3s_N[0, 1], I_3_3s_N[1, 2], I_3_3s_N[2, 0])

   I_1_O = I_1_1s + I_1s_O
   I_2_O = I_2_2s + I_2s_O
   I_3_O = I_3_3s + I_3s_O
   print("\nI_1_rel_O = {0}".format(I_1_O))
   print("\nI_2_rel_O = {0}".format(I_2_O))
   print("\nI_3_rel_O = {0}".format(I_3_O))

   # assembly inertia tensor is the sum of the inertia tensor of each component
   I_B_O = I_1_O + I_2_O + I_3_O
   print("\nI_B_rel_O = {0}".format(I_B_O))

   # n_a is parallel to line L
   n_a = sqrt(2) / 2 * (N.x + N.z)
   print("\nn_a = {0}".format(n_a))

   # calculate moment of inertia of for point O of assembly about line L
   I_a_a_B_O = simplify(dot(n_a, dot(I_B_O, n_a)))
   print("\nI_a_a_B_rel_O = {0}".format(I_a_a_B_O))

   # use the second value since k_a is non-negative
   k_a_val = solve(I_a_a_B_O - 3 * m * k_a**2, k_a)[1]
   print("\nk_a = {0}".format(k_a_val))

   dk_a_dtheta = diff(k_a_val, theta)
   print("\ndk_a/dtheta = {0}".format(dk_a_dtheta))

   # solve for theta = 0 using a simplified expression or
   # else no solution will be found
   soln = solve(3*cos(theta) + 12*sin(theta) - 4*sin(theta)*cos(theta), theta)
   # ignore complex values of theta
   theta_vals = [s for s in soln if s.is_real]

   print("")
   for th in theta_vals:
       print("k_a({0}) = {1}".format((th * 180 / pi).evalf(),
                                     k_a_val.subs(theta, th).evalf()))

Exercise 6.10
=============

.. jupyter-kernel:: python3
   :id: ex6-10

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex6-10` or Jupyter notebook:
   :jupyter-download:notebook:`ex6-10`.

.. jupyter-execute::

   from sympy import Matrix
   from sympy import pi, acos
   from sympy import simplify, symbols
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import inertia, inertia_of_point_mass
   from sympy.physics.mechanics import dot


   def inertia_matrix(dyadic, rf):
      """Return the inertia matrix of a given dyadic for a specified
      reference frame.
      """
      return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])


   def angle_between_vectors(a, b):
      """Return the minimum angle between two vectors. The angle returned for
      vectors a and -a is 0.
      """
      angle = (acos(dot(a, b)/(a.magnitude() * b.magnitude())) *
               180 / pi).evalf()
      return min(angle, 180 - angle)


   m, m_R, m_C, rho, r = symbols('m m_R m_C rho r', real=True, nonnegative=True)

   N = ReferenceFrame('N')
   pA = Point('A')
   pPs = pA.locatenew('P*', 3*r*N.x - 2*r*N.y)

   m_R = rho * 24 * r**2
   m_C = rho * pi * r**2
   m = m_R - m_C

   I_Cs_A = inertia_of_point_mass(m, pPs.pos_from(pA), N)
   I_C_Cs = inertia(N, m_R*(4*r)**2/12 - m_C*r**2/4,
                     m_R*(6*r)**2/12 - m_C*r**2/4,
                     m_R*((4*r)**2+(6*r)**2)/12 - m_C*r**2/2)

   I_C_A = I_C_Cs + I_Cs_A
   print("\nI_C_rel_A = {0}".format(I_C_A))

   # Eigenvectors of I_C_A are the parallel to the principal axis for point A
   # of Body C.
   evecs_m = [triple[2]
             for triple in inertia_matrix(I_C_A, N).eigenvects()]

   # Convert eigenvectors from Matrix type to Vector type.
   evecs = [sum(simplify(v[0][i]).evalf() * n for i, n in enumerate(N))
            for v in evecs_m]

   # N.x is parallel to line AB
   print("\nVectors parallel to the principal axis for point A of Body C and the" +
         "\ncorresponding angle between the principal axis and line AB (degrees):")
   for v in evecs:
       print("{0}\t{1}".format(v, angle_between_vectors(N.x, v)))

Exercise 6.13
=============

.. jupyter-kernel:: python3
   :id: ex6-13

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex6-13` or Jupyter notebook:
   :jupyter-download:notebook:`ex6-13`.

.. jupyter-execute::

   from sympy import S, Matrix
   from sympy import pi, acos
   from sympy import simplify, sqrt, symbols
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import inertia_of_point_mass
   from sympy.physics.mechanics import dot
   import numpy as np


   def inertia_matrix(dyadic, rf):
      """Return the inertia matrix of a given dyadic for a specified
      reference frame.
      """
      return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])


   def convert_eigenvectors_matrix_vector(eigenvectors, rf):
      """Return a list of Vectors converted from a list of Matrices.
      rf is the implicit ReferenceFrame for the Matrix representation of the
      eigenvectors.
      """
      return [sum(simplify(v[0][i]).evalf() * n for i, n in enumerate(N))
               for v in eigenvectors]


   def angle_between_vectors(a, b):
      """Return the minimum angle between two vectors. The angle returned for
      vectors a and -a is 0.
      """
      angle = (acos(dot(a, b) / (a.magnitude() * b.magnitude())) *
               180 / pi).evalf()
      return min(angle, 180 - angle)


   m = symbols('m', real=True, nonnegative=True)
   m_val = 1
   N = ReferenceFrame('N')
   pO = Point('O')
   pP = pO.locatenew('P', -3 * N.y)
   pQ = pO.locatenew('Q', -4 * N.z)
   pR = pO.locatenew('R', 2 * N.x)
   points = [pO, pP, pQ, pR]

   # center of mass of assembly
   pCs = pO.locatenew('C*', sum(p.pos_from(pO) for p in points) / S(len(points)))
   print(pCs.pos_from(pO))

   I_C_Cs = (inertia_of_point_mass(m, points[0].pos_from(pCs), N) +
            inertia_of_point_mass(m, points[1].pos_from(pCs), N) +
            inertia_of_point_mass(m, points[2].pos_from(pCs), N) +
            inertia_of_point_mass(m, points[3].pos_from(pCs), N))
   print("I_C_Cs = {0}".format(I_C_Cs))

   # calculate the principal moments of inertia and the principal axes
   M = inertia_matrix(I_C_Cs, N)

   # use numpy to find eigenvalues/eigenvectors since sympy failed
   # note that the eigenvlaues/eigenvectors are the
   # prinicpal moments of inertia/principal axes
   eigvals, eigvecs_np = np.linalg.eigh(np.matrix(M.subs(m, m_val).n().tolist(), dtype=float))
   eigvecs = [sum(eigvecs_np[i, j] * n for i, n in enumerate(N))
            for j in range(3)]

   # get the minimum moment of inertia and its associated principal axis
   e, v = min(zip(eigvals, eigvecs))

   # I = m * k**2, where I is the moment of inertia,
   # m is the mass of the body, k is the radius of gyration
   k = sqrt(e / (4 * m_val))
   print("\nradius of gyration, k = {0} m".format(k))

   # calculate the angle between the associated principal axis and the line OP
   # line OP is parallel to N.y
   theta = angle_between_vectors(N.y, v)
   print("\nangle between associated principal axis and line OP = {0}°".format(theta))

Exercise 6.14
=============

.. jupyter-kernel:: python3
   :id: ex6-14

.. note::
   You can download this example as a Python script:
   :jupyter-download:script:`ex6-14` or Jupyter notebook:
   :jupyter-download:notebook:`ex6-14`.

.. jupyter-execute::

   from sympy import Matrix, S
   from sympy import integrate, pi, simplify, solve, symbols
   from sympy.physics.mechanics import ReferenceFrame, Point
   from sympy.physics.mechanics import cross, dot
   from sympy.physics.mechanics import inertia, inertia_of_point_mass


   def inertia_matrix(dyadic, rf):
      """Return the inertia matrix of a given dyadic for a specified
      reference frame.
      """
      return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])


   def integrate_v(integrand, rf, bounds):
      """Return the integral for a Vector integrand."""
      return sum(simplify(integrate(dot(integrand, n), bounds)) * n for n in rf)


   def index_min(values):
      return min(range(len(values)), key=values.__getitem__)


   m, a, b, c = symbols('m a b c', real=True, nonnegative=True)
   x, y, r = symbols('x y r', real=True)
   k, n = symbols('k n', real=True, positive=True)
   N = ReferenceFrame('N')

   # calculate right triangle density
   V = b * c / 2
   rho = m / V

   # Kane 1985 lists scalars of inertia I_11, I_12, I_22 for a right triangular
   # plate, but not I_3x.
   pO = Point('O')
   pCs = pO.locatenew('C*', b/3*N.x + c/3*N.y)
   pP = pO.locatenew('P', x*N.x + y*N.y)
   p = pP.pos_from(pCs)

   I_3 = rho * integrate_v(integrate_v(cross(p, cross(N.z, p)),
                                    N, (x, 0, b*(1 - y/c))),
                        N, (y, 0, c))
   print("I_3 = {0}".format(I_3))

   # inertia for a right triangle given ReferenceFrame, height b, base c, mass
   inertia_rt = lambda rf, b_, c_, m_: inertia(rf,
                  m_*c_**2/18,
                  m_*b_**2/18,
                  dot(I_3, N.z),
                  m_*b_*c_/36,
                  dot(I_3, N.y),
                  dot(I_3, N.x)).subs({m:m_, b:b_, c:c_})

   theta = (30 + 90) * pi / 180
   N2 = N.orientnew('N2', 'Axis', [theta, N.x])

   # calculate the mass center of the assembly
   # Point O is located at the right angle corner of triangle 1.
   pCs_1 = pO.locatenew('C*_1', 1/S(3) * (a*N.y + b*N.x))
   pCs_2 = pO.locatenew('C*_2', 1/S(3) * (a*N2.y + b*N2.x))
   pBs = pO.locatenew('B*', 1/(2*m) * m * (pCs_1.pos_from(pO) +
                                          pCs_2.pos_from(pO)))
   print("\nB* = {0}".format(pBs.pos_from(pO).express(N)))

   # central inertia of each triangle
   I1_C_Cs = inertia_rt(N, b, a, m)
   I2_C_Cs = inertia_rt(N2, b, a, m)

   # inertia of mass center of each triangle about Point B*
   I1_Cs_Bs = inertia_of_point_mass(m, pCs_1.pos_from(pBs), N)
   I2_Cs_Bs = inertia_of_point_mass(m, pCs_2.pos_from(pBs), N)

   I_B_Bs = I1_C_Cs + I1_Cs_Bs + I2_C_Cs + I2_Cs_Bs
   print("\nI_B_Bs = {0}".format(I_B_Bs.express(N)))

   # central principal moments of inertia
   evals = inertia_matrix(I_B_Bs, N).eigenvals()

   print("\neigenvalues:")
   for e in evals.keys():
      print(e)

   print("\nuse a/b = r")
   evals_sub_a = [simplify(e.subs(a, r*b)) for e in evals.keys()]
   for e in evals_sub_a:
       print(e)

   for r_val in [2, 1/S(2)]:
       print("\nfor r = {0}".format(r_val))
       evals_sub_r = [e.subs(r, r_val) for e in evals_sub_a]
       print("eigenvalues:")
       for e in evals_sub_r:
           print("{0} = {1}".format(e, e.n()))

       # substitute dummy values for m, b so that min will actually work
       min_index = index_min([e.subs({m : 1, b : 1}) for e in evals_sub_r])
       min_e = evals_sub_r[min_index]
       print("min: {0}".format(min_e))

       k_val = solve(min_e - 2*m*k**2, k)
       assert(len(k_val) == 1)
       print("k = {0}".format(k_val[0]))
       n_val = solve(k_val[0] - n*b, n)
       assert(len(n_val) == 1)
       print("n = {0}".format(n_val[0]))

   print("\nResults in text: n = 1/3, sqrt(35 - sqrt(241))/24")
