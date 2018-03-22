PyDy for Autolev Users
==========================

Introduction
----------------

Autolev (now defunct) is a domain specific programming language which is used for
symbolic multibody dynamics. The SymPy mechanics module now has enough
power and functionality to be a fully featured symbolic dynamics module.
The PyDy package extends the SymPy output to the numerical domain for
simulation, analyses and visualization. Autolev and PyDy have a lot in
common but there are also many differences between them. This page shall
expand upon their differences. It is meant to be a go-to reference for
Autolev users who want to transition to PyDy.

It would be helpful to have a basic idea about scientific computing with
Python before going over this page. If you aren’t familiar with
scientific computing with Python there are many sources to learn. You
can start with the Python programming language itself, with the
canonical source being the `Python
Documentation <https://docs.python.org/>`__. The `SciPy
Lectures <http://www.scipy-lectures.org/>`__ are a great intro to
scientific computing with Python. Finally, to learn about how to do
symbolic manipulation with SymPy, check out the `SymPy
Documentation <http://docs.sympy.org/>`__, especially the tutorial.

Some Key Differences
------------------------

+-----------------------------------+-----------------------------------+
|          **Autolev**              |             **PyDy**              |            
+===================================+===================================+
||                                  ||                                  | 
| Autolev is a domain specific      | PyDy is a library written in the  |
| programming language designed to  | general purpose language Python.  |
| perform multibody dynamics. Since | Although Autolev's code is more   |
| it is a language of its own, it   | compact, PyDy(by virtue of being  |
| has a very rigid language         | an add on to Python) is more      |
| specification. It predefines and  | flexible. The users have more     |
| computes many things based on the | control over what they can do. For|
| input code. Its statements are a  | example, one can create a class in|
| lot cleaner as a result of this.  | their code for let's say a type of|
|                                   | rigibodies with with common       |
|                                   | properties.                       |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev generates Matlab, C, or   | PyDy generates numerical Python,  |
| Fortan code from a small set of   | C or Octave/Matlab code from a    |
| symbolic mathematics.             | large set of symbolic mathematics |
|                                   | created with SymPy. It also builds|
|                                   | on the popular scientific Python  |
|                                   | stack such as NumPy, SciPy,       |
|                                   | IPython, matplotlib, Cython and   |
|                                   | Theano.                           |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev uses 1 (one) based        | Python uses 0 (zero) based        |
| indexing. The initial element of  | indexing. The initial element of  |
| a sequence is found using a[1].   | a sequence is found using a[0].   |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev is case insensitive.      | PyDy code being Python code is    |
|                                   | case sensitive.                   |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| One can define their own commands | PyDy code is Python code, so one  |
| in Autolev by making .R and .A    | can define functions in their     |
| files which can be used in their  | code. This is a lot more          |
| programs.                         | convenient.                       |
+-----------------------------------+-----------------------------------+
||                                  ||                                  |
| Autolev is proprietary.           | PyDy is open source.              |
+-----------------------------------+-----------------------------------+

Rough Autolev-PyDy Equivalents
----------------------------------

The tables below give rough equivalents for some common Autolev
expressions. **These are not exact equivalents**, but rather should be
taken as hints to get you going in the right direction. For more detail
read the built-in documentation on `SymPy vectors <http://docs.sympy.org/latest/modules/physics/vector/index.html>`__
, `SymPy mechanics <http://docs.sympy.org/latest/modules/physics/mechanics/index.html>`__ and
`PyDy <http://www.pydy.org/documentation.html>`__ .

In the tables below, it is assumed that you have executed the following
commands in Python:
::

	import sympy.physics.mechanics as me
	import sympy as sm

Mathematical Equivalents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-----------------------+-----------------------+-----------------------+
| **Autolev**           | **PyDy**              | **Notes**             |
+=======================+=======================+=======================+
||                      ||                      ||                      |
| Constants A, B        | a, b = sm.symbols     | Note that the names   |
|                       | (‘a b’, real=True)    | of the symbols can be |
|                       |                       | different from the    |
|                       |                       | names of the          |
|                       |                       | variables they are    |
|                       |                       | assigned to. We can   |
|                       |                       | define a, b =         |
|                       |                       | symbols(‘b a’) but it |
|                       |                       | is good practice to   |
|                       |                       | follow the            |
|                       |                       | convention.           |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Constants C+          | c = sm.symbols(‘c’,   | Refer to SymPy        |
|                       | real=True, nonnegative| assumptions for more  |
|                       | =True)                | information.          |
|                       |                       | `modules/core         |
|                       |                       | <http://docs.sy       |
|                       |                       | mpy.org/latest/module |
|                       |                       | s/core.html>`__       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Constants D-          | d = sm.symbols(‘d’,   |                       |
|                       | real=True,            |                       |
|                       | nonpositive=True)     |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Specified Phi         | Phi =                 |                       | 
|                       | me.dynamicsymbols(‘Phi|                       |
|                       | ')                    |                       |       
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Variables q, s        | q, s =                |                       |
|                       | me.dynamicsymbols     |                       |
|                       | (q, s)                |                       |
+-----------------------+-----------------------+-----------------------+
| Variables x’’         | x =                   |                       |
|                       | me.dynamicsymbols(‘x’)|                       |
|                       |                       |                       |
|                       | xd =                  |                       |
|                       | me.dynamicsymbols     |                       |
|                       | (‘x’, 1)              |                       |
|                       |                       |                       |
|                       | xd2 =                 |                       |
|                       | me.dynamicsymbols(‘x’,|                       |
|                       | 2)                    |                       |
+-----------------------+-----------------------+-----------------------+
| Variables y{2}’       | y1 =                  |                       |
|                       | me.dynamicsymbols     |                       |
|                       | (‘y1’)                |                       |
|                       |                       |                       |
|                       | y2 =                  |                       |
|                       | me.dynamicsymbols     |                       |
|                       | (‘y2’)                |                       |
|                       |                       |                       |
|                       | y1d =                 |                       |
|                       | me.dynamicsymbols(‘y1’|                       |
|                       | , 1)                  |                       |
|                       |                       |                       |
|                       | y2d =                 |                       |
|                       | me.dynamicsymbols(‘y2 |                       |
|                       | , 1)                  |                       |
+-----------------------+-----------------------+-----------------------+
| MotionVariables u{2}  |u1 =                   | SymPy doesn’t         |
|                       |me.dynamicsymbols(‘u1’)| differentiate between |
|                       |                       | variables,            |
|                       |u2 =                   | motionvariables and   |
|                       |me.dynamicsymbols(‘u2’)| specifieds during     |
|                       |                       | declaration. Instead, |
|                       |                       | it takes different    |
|                       |                       | lists of these as     |
|                       |                       | parameters in objects |
|                       |                       | like the KanesMethod. |
+-----------------------+-----------------------+-----------------------+
| Imaginary j           | j = sm.symbols(‘j’)   | I is a sympy object   |
|                       |                       | which stands for the  |
|                       | j = I                 | imaginary unit. One   |
|                       |                       | can define complex    |
|                       |                       | numbers using it.     |
|                       |                       |                       |
|                       |                       | z = x + I*y           |
|                       |                       |                       |
|                       |                       | where x, y and z are  |
|                       |                       | symbols.              |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Tina = 2*pi           | Tina = 2*sm.pi        | If one wants to use   |
|                       |                       | numerical constants   |
|                       |                       | instead of symbolic   |
|                       |                       | ones they can import  |
|                       |                       | them from mpmath.     |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| abs(x)^3 + sin(x)^2 + | sm.abs(x)**3          |                       |
| acos(x)               | + sm.sin(x)**2+       |                       |
|                       | + sm.acos(x)          |                       |
+-----------------------+-----------------------+-----------------------+
| E = (x+2*y)^2 +       | E = (x+2*y)**2 +      | For more information  |
| 3*(7+x)*(x+y)         | 3*(7+x)*(x+y)         | refer to              | 
|                       |                       | `simplification. <htt |
| Expand(E, n:m)        | sm.expand(E)          | p://docs.sympy.org/la |
|                       |                       | test/tutorial/simplif |
| Factor(E, x)          | sm.factor(E)          | ication.html>`__      |
|                       |                       |                       |
| Coef(y, x)            | y.coeff(x)            |                       |
|                       |                       |                       |
| Replace(y, sin(x)=3)  | y = y.xreplace        |                       |
|                       | ({sm.sin(x): 3})      |                       |
+-----------------------+-----------------------+-----------------------+
| Dy = D(E, y)          | sm.diff(E, y)         | For more information  |
|                       |                       | refer to `calculus.   |
| Dt = Dt(E)            |                       | <http: //docs.sympy.or|
|                       |                       | g/latest/tutorial/    |
|                       | sm.diff(E, t) where t | calculus.html>`__     |
|                       | = me.dynamicsymbols._t|                       |
|                       |                       |                       |
|                       |                       |                       |
|                       | Works if the          |                       |
|                       | expression is made up |                       |
|                       | of dynamicsymbols.    |                       |
+-----------------------+-----------------------+-----------------------+
| TY = Taylor(x*cos(x), | ty =                  | Execute               |
| 0:7, x = 0)           | taylor(x*sm.cos(x)    |                       |
|                       | ,0 , 7)               | from                  |
|                       |                       | sympy.mpmath.import \*|             
|                       |                       |                       |
|                       |                       | For more information  |
|                       |                       | refer to              |
|                       |                       | `mpmath/calculus. <ht |
|                       |                       | tp://docs             |
|                       |                       | .sympy.org/0.7.6/modu |
|                       |                       | les/mpmath/calculus/a |
|                       |                       | pproximation.html>`__ |
+-----------------------+-----------------------+-----------------------+
| F = Evaluate(E, x=a,  | E.subs([(x, a), (y,   |                       |
| y=2)                  | 2)])                  |                       |
|                       |                       |                       |
|                       | To get floating point |                       |
|                       | numbers from numerical|                       |
|                       | expressions use evalf |                       |
|                       |                       |                       |
|                       | E.evalf((a+sm.pi).subs|                       |
|                       | ({a:3}))              |                       |                  
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| P = Polynomial([a, b, | p = sm.Poly(a*x**2    | For more information  |
| c], x)                | + b*x + c)            | refer to              |
|                       |                       | `modules/polys. <htt  |
|                       |                       | p://docs.sympy.org/la |
|                       |                       | test/modules/polys/re |
|                       |                       | ference.html>`__      |
+-----------------------+-----------------------+-----------------------+
| Roots(Polynomial([a,  | sm.solve(             | For more information  |
| b, c], x), x, 2)      | sm.Poly(a*x**2 +      | refer to `solvers. <ht| 
|                       | b*x + c))             | tp://docs.sympy.org/la|
|                       |                       | test/modules/solvers/ |
|                       |                       | solvers.html>`__      |
|                       |                       |                       |
|                       |                       | For numerical         |
|                       |                       | computation related   |
|                       |                       | to polynomials and    |
|                       |                       | roots refer to        |
|                       |                       | `mpmath/calculus. <htt|
|                       |                       | p://docs.s            | 
|                       |                       | ympy.org/0.7.6/module |
|                       |                       | s/mpmath/calculus/pol |
|                       |                       | ynomials.html>`__     |
+-----------------------+-----------------------+-----------------------+
| Solve(A, x1, x2)      | sm.linsolve(A,        | For more information  |
|                       | (x1, x2))             | refer to              |   
|                       |                       | `solvers/solveset. <ht|
| where A is an         | where A is an         | tp://docs.sympy.org/l |
| augmented matrix that | augmented matrix      | atest/modules/solvers |
| represents the linear |                       | /solveset.html>`__    |
| equations and x1, x2  |                       |                       |
| are the variables to  |                       |                       |
| solve for.            |                       |                       |
|                       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| RowMatrix = [1, 2, 3, | row_matrix =          | For more information  |
| 4]                    | Matrix([[1],[2],      | refer to `matrices. <h|
|                       | [3],[4]])             | ttp://docs.sympy.org/ |
|                       |                       | latest/tutorial/      |            
|                       | col_matrix =          | matrices.html>`__     |                     
| ColMatrix = [1; 2; 3; | Matrix([1, 2, 3, 4])  |                       |           
| 4]                    |                       |                       |
|                       | MO = Matrix([[a, b],  |                       |
| MO = [a, b; c, 0]     | [c, 0]])              |                       |
|                       |                       |                       |
| MO[2, 2] := d         | MO[1, 1] = d          |                       |
|                       |                       |                       |
| A + B*C               | A + B*C               |                       |
|                       |                       |                       |
|                       | A.shape(0)            |                       |
|                       |                       |                       |
|                       | A.shape(1)            |                       |
| Cols(A)               |                       |                       |
|                       | M.det()               |                       |
| Rows(A)               |                       |                       |
|                       | M[2, 3]               |                       |
| Det(A)                |                       |                       |
|                       | M**-1                 |                       |
| Element(A, 2, 3)      |                       |                       |
|                       | trace(A)              |                       |
| Inv(A)                |                       |                       |
|                       | A.T                   |                       |
| Trace(A)              |                       |                       |
|                       | A.eigenvals()         |                       |
| Transpose(A)          |                       |                       |
|                       |                       |                       |
| Eig(A)                |                       |                       |
+-----------------------+-----------------------+-----------------------+


Physical Equivalents
~~~~~~~~~~~~~~~~~~~~~~~~

+-----------------------+-----------------------+-----------------------+
| **Autolev**           | **PyDy**              | **Notes**             |
+=======================+=======================+=======================+
| Bodies A              | m = symbol(‘m’)       | The 4th and 5th       |
|                       |                       | arguments are for the |
| Declares A, its       | Ao = symbols(‘Ao’)    | mass and inertia.     |
| masscenter Ao, and    |                       | These are specified   |
| orthonormal vectors   | Af =                  | after the declaration |
| A1>, A2> and A3>      | ReferenceFrame(‘Af’)  | in Autolev.           |
| fixed in A.           |                       |                       |
|                       | I = outer(Af.x, Af.x) | One can pass in None  |
|                       |                       | for the parameters    |
|                       | P = Point(‘P’)        | and use setters       |
|                       |                       | A.mass = \_ and       |
|                       | A = RigidBody(‘A’,    | A.inertia = \_ to set |
|                       | Ao, Af, m, (I, P))    | them later.           |
|                       |                       |                       |
|                       | Af.x, Af.y and Af.z   | For more information  |
|                       | are equivalent to     | refer to              |
|                       | A1>, A2> and A3>.     | `mechanics/masses. <ht|
|                       |                       | tp://docs.sym         |
|                       |                       | py.org/latest/modules |
|                       |                       | /physics/mechanics/ma |
|                       |                       | sses.html>`__         |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Frames B              | B =                   | For more information  |
|                       | ReferenceFrame(‘B’)   | refer to              |
|                       |                       | `physics/vector. <http|
|                       |                       | ://docs.sympy         |
|                       |                       | .org/latest/modules/p |
|                       |                       | hysics/vector/vectors |
|                       |                       | .html>`__             |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Newtonian N           | N =                   | SymPy doesn’t specify |
|                       | ReferenceFrame(‘N’)   | that a frame is       |
|                       |                       | inertial during       |
|                       |                       | declaration. Many     |
|                       |                       | functions such as     |
|                       |                       | set_ang_vel() take    |
|                       |                       | the inertial          |
|                       |                       | reference frame as a  |
|                       |                       | parameter.            |
+-----------------------+-----------------------+-----------------------+
| Particles C           | m = symbol(‘m’)       | The 2nd and 3rd       |
|                       |                       | arguments are for the |
|                       | po = Point(‘po’)      | point and mass. In    |
|                       |                       | Autolev, these are    |
|                       | C = Particle(‘C’, po, | specified after the   |
|                       | m)                    | declaration..         |
|                       |                       |                       |
|                       |                       | One can pass in None  |
|                       |                       | and use setters       |
|                       |                       | (A.point = \_ and     |
|                       |                       | A.mass = \_) to set   |
|                       |                       | them later.           |
+-----------------------+-----------------------+-----------------------+
| Points P, Q           | P = Point(‘P’)        |                       |
|                       |                       |                       |
|                       | Q = Point(‘Q’)        |                       |
+-----------------------+-----------------------+-----------------------+
| Mass B=mB             | mB = symbols(‘mB’)    |                       |
|                       |                       |                       |
|                       | B.mass = mB           |                       |
+-----------------------+-----------------------+-----------------------+
| Inertia B,I1,I2,I3,I12|I = inertia(Bf, i1, i2,| For more information  |
| I23,I31               |i3, i12, i23, i31)     | refer to the          |
|                       |                       | `mechanics api. <http:|
|                       |B.inertia = (I, P)     | //docs.sympy.org/lates|
|                       |where B is a rigidbody,| t/modules/physics/mech|
|                       |Bf is the related frame| anics/api/part_bod.   |
|                       |and P is the center of | html>`__              |
|                       |mass of B.             |                       |
|                       |                       |                       |
|                       |Inertia dyadics can    |                       |
|                       |also be formed using   |                       |
|                       |vector outer products. |                       |
|                       |                       |                       |
|                       |I = outer(N.x, N.x)    |                       |
+-----------------------+-----------------------+-----------------------+
| vec> = P_O_Q>/L       | vec  = (Qo.pos_from   |For more information   |                
|                       | (O))/L                |refer to               |
| vec> = u1*N1> + u2*N2>| vec = u1*N.x + u2*N.y |`physics/vector. <http:|
|                       |                       |//docs.sympy.org/latest|
| Cross(a>, b>)         | cross(a, b) where a   |modules/physics/vector |
|                       | and b are vectors     |/index.html>`__        |
|                       |                       |                       |
| Dot(a>, b>)           | dot(a, b)             |                       |
|                       |                       |                       |
| Mag(v>)               | v.magnitude()         |                       |
|                       |                       |                       |
| Unitvec(v>)           | v.normalize()         |                       |
+-----------------------+-----------------------+-----------------------+
| P_O_Q> = LA*A1>       | Q.point =             | For more information  |
|                       | O.locatenew(‘Qo’,     | refer to the          | 
| where O is a point    | LA*A.x)               | `kinematics api. <http|
|                       |                       | ://docs.sympy.org/late|
| P_P_Q> = LA*A1>       | where A is a          | st/modules/physics/vec|
|                       | reference frame.      | tor/api/kinematics.   |
| where P is a particle |                       | html>`__              |
|                       | Q.point =             |                       |   
|                       | P.point.locatenew(‘Qo |                       | 
|                       | ’,                    |                       |
|                       | LA*A.x)               |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| V_O_N> = u3*N.1> +    | O.set_vel(N, u1*N.x + |                       |
| u4*N.2>               | u2*N.y)               |                       |
+-----------------------+-----------------------+-----------------------+
| A_O_N> = 0>           | O.set_acc(N, 0)       |                       |
|                       |                       |                       |
| Acceleration of point |                       |                       |
| O in reference frame  |                       |                       |
| N.                    |                       |                       |
+-----------------------+-----------------------+-----------------------+
| W_B_N> = qB’*B3>      | B.set_ang_vel(N,      |                       |
|                       | qBd*Bf.z)             |                       |
| Angular velocity of   |                       |                       |
| body B in reference   | where Bf is the frame |                       |
| frame F.              | associated with the   |                       |
|                       | body B.               |                       |
+-----------------------+-----------------------+-----------------------+
| ALF_B_N> = Dt(W_B_N>, | B.set_ang_acc(N,      |                       |
| N)                    | diff(B.ang_vel_in(N)) |                       |
|                       |                       |                       |
| Angular acceleration  |                       |                       |
| of body B in          |                       |                       |
| reference frame N.    |                       |                       |
+-----------------------+-----------------------+-----------------------+
| Force_O> = F1*N1> +   | In SymPy one should   |                       |
| F2*N2>                | have a list which     |                       |
|                       | contains all the      |                       |
| Torque_A> =           | forces and torques.   |                       |
| -c*qA’*A3>            |                       |                       |
|                       | fL.append((O, f1*N.x  |                       |
|                       | + f2*N.y))            |                       |
|                       |                       |                       |
|                       | where fL is the force |                       |
|                       | list.                 |                       |
|                       |                       |                       |
|                       | fl.append((A,         |                       |
|                       | -c*qAd*A.z))          |                       |
+-----------------------+-----------------------+-----------------------+
| A_B                   | A.dcm(B)              |                       |
|                       |                       |                       |
| or                    |                       |                       |
|                       |                       |                       |
| Dircos(A,B)           |                       |                       |
|                       |                       |                       |
| where A and B are     |                       |                       |
| reference frames      |                       |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| CM(B)                 | B.masscenter          |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Mass(A,B,C)           | A.mass + B.mass +     |                       |
|                       | C.mass                |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| V1pt(A,B,Bq,Q)        | Q.v1pt_theory(Bq, A,  |                       |
|                       | B)                    |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| V2pts(A,B,P,Q)        | Q.v2pt_theory(P, A,   |                       |
|                       | B)                    |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| A1pt(A,B,Bq,Q)        | Q.a1pt_theory(Bq, A,  |                       |
|                       | B)                    |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| A2pts(A,B,P,Q)        | Q.a2pt_theory(P, A,   |                       |
|                       | B)                    |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Angvel(A,B)           | B.ang_vel_in(A)       |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Simprot(A, B, 1, x)   | B.orient(A, ‘Axis’,   |                       |
|                       | x, A.x)               |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Gravity(G*N1>)        | fL.extend(gravity(g*N |                       |
|                       | .x,                   |                       |
|                       | P1, P2, ...))         |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Force(P/Q, v>)        | fL.append((P, -1*v),  |                       |
|                       | (Q, v))               |                       |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Torque(A/B, v>)       | fL.append((A, -1*v),  |                       |
|                       | (B, v))               |                       |
+-----------------------+-----------------------+-----------------------+
| Fr()                  | (fr, frstar) =        | For more details      |
|                       | KM.kanes_equations(fL | refer to              |
| FrStar()              | ,                     | `mechanics/kane <http:|
|                       | bL)                   | //docs.sympy.org/lates|
|                       |                       | t/modules/physics/mech|
|                       | where KM is the       | anics/kane.html>`__   |
|                       | KanesMethod object.   | and                   |                       
|                       |                       | `the api. <http://docs| 
|                       |                       | .sympy.org/0.7.5/modul|
|                       |                       | es/physics/mechanics/a|
|                       |                       | pi/kane.html>`__      |
+-----------------------+-----------------------+-----------------------+
||                      ||                      ||                      |
| Kindiffs(A, B ...)    | KM.kindiffdict()      |                       |
+-----------------------+-----------------------+-----------------------+
| Momentum(option)      | linear_momentum(N,    |                       |
|                       | B1, B2 ...)           |                       |
|                       |                       |                       |
|                       | reference frame       |                       |
|                       | followed by one or    |                       |
|                       | more bodies           |                       |
|                       |                       |                       |
|                       | angular_momentum(O,   |                       |
|                       | N, B1, B2 ...)        |                       |
|                       |                       |                       |
|                       | point, reference      |                       |
|                       | frame followed by one |                       |
|                       | or more bodies        |                       |
+-----------------------+-----------------------+-----------------------+
| KE()                  | kinetic_energy(N, B1, |                       |
|                       | B2 ...)               |                       |
|                       |                       |                       |
|                       | reference frame       |                       |
|                       | followed by one or    |                       |
|                       | more bodies           |                       |
+-----------------------+-----------------------+-----------------------+
| Constrain(...)        | velocity_constraints  | For more details      |
|                       | = [...]               | refer to              |
|                       |                       | `mechanics/kane <http |
|                       | u_dependent = [...]   | ://docs.sympy.or      |
|                       |                       | g/latest/modules/phys |
|                       | u_auxiliary = [...]   | ics/mechanics/kane.ht |
|                       |                       | ml>`__ and the        |
|                       | These lists are       | `kane api. <htt       |
|                       | passed to the         | p://docs.sympy.org/0.7|
|                       | KanesMethod object    | .5/modules/physics/mec|
|                       |                       | hanics/api/kane.      |
|                       |                       | html>`__              |
|                       |                       |                       |
|                       |                       |                       |
|                       |                       |                       |
+-----------------------+-----------------------+-----------------------+
| Kane()                |KanesMethod(frame,     | For more details      |
|                       |q_ind,u_ind,kd_eqs,    | refer to              |
|                       |q_dependent,configura  | `mechanics/kane <http |
|                       |tion_constraints,u_dep | ://docs.sympy.or      |
|                       |endent,velocity_const  | g/latest/modules/phys |
|                       |raints,acceleration_c  | ics/mechanics/kane.ht |
|                       |onstraints,u_auxiliary)| ml>`__ and the        |
|                       |                       | `kane api. <htt       |
|                       |The KanesMethod        | p://docs.sympy.org/0.7|
|                       |object takes a         | .5/modules/physics/mec|
|                       |reference frame        | hanics/api/kane.      |
|                       |followed by multiple   | html>`__              |
|                       |lists as arguments.    |                       |
|                       |                       |                       |
|                       |                       |                       |
|                       |                       |                       |
|                       |                       |                       |
|                       |                       |                       |
+-----------------------+-----------------------+-----------------------+

Numerical Evaluation and Visualization
----------------------------------------

Autolev’s CODE Option() command allows one to generate Matlab, C, or
Fortran code for numerical evaluation and visualization. Option can be
Dynamics, ODE, Nonlinear or Algebraic.

Numerical evaluation for dynamics can be achieved using PyDy. One can
pass in the KanesMethod object to the System class along with the values
for the constants, specifieds, initial conditions and time steps. The
equations of motion can then be integrated. The plotting is achieved
using matlplotlib. Here is an example from the PyDy documentation on how
it is done:
::
	  from numpy import array, linspace, sin
	  from pydy.system import System

	  sys = System(kane,
	  		    constants = {mass: 1.0, stiffness: 1.0,
	  		                 damping: 0.2, gravity: 9.8},
	  		    specifieds = {force: lambda x, t: sin(t)},
	  		    initial_conditions = {position: 0.1, speed:-1.0},
	  		    times = linspace(0.0, 10.0, 1000))

	  y = sys.integrate()

	  import matplotlib.pyplot as plt
	  plt.plot(sys.times, y)
	  plt.legend((str(position), str(speed)))
	  plt.show()

For information on all the things PyDy can accomplish refer to the PyDy
documentation.

The tools in the PyDy workflow are :

-  `SymPy <http://sympy.org>`__: SymPy is a Python library for
    symbolic computation. It provides computer algebra capabilities
    either as a standalone application, as a library to other
    applications, or live on the web as SymPy Live or SymPy Gamma.

-  `NumPy <http://www.numpy.org/>`__: NumPy is a library for the
    Python programming language, adding support for large,
    multi-dimensional arrays and matrices, along with a large
    collection of high-level mathematical functions to operate on
    these arrays.

-  `SciPy <https://www.scipy.org/>`__: SciPy is an open source
    Python library used for scientific computing and technical
    computing. SciPy contains modules for optimization, linear
    algebra, integration, interpolation, special functions, FFT,
    signal and image processing, ODE solvers and other tasks common
    in science and engineering.

-  `IPython <https://ipython.org/>`__: IPython is a command shell
    for interactive computing in multiple programming languages,
    originally developed for the Python programming language, that
    offers introspection, rich media, shell syntax, tab completion,
    and history.

-  `Theano <http://deeplearning.net/software/theano/>`__: Theano is
    a numerical computation library for Python. In Theano,
    computations are expressed using a NumPy-esque syntax and
    compiled to run efficiently on either CPU or GPU architectures

-  `Cython <http://cython.org/>`__: Cython is a superset of the
    Python programming language, designed to give C-like performance
    with code that is mostly written in Python. Cython is a compiled
    language that generates CPython extension modules.

-  `matplotlib <https://matplotlib.org/>`__: matplotlib is a
    plotting library for the Python programming language and its
    numerical mathematics extension NumPy.

One will be able to write code equivalent to the Matlab, C or Fortran
code generated by Autolev using these scientific computing tools. It is
recommended to go over these modules to gain an understanding of
scientific computing with Python.

Links
----------

`SymPy tutorial <http://docs.sympy.org/latest/tutorial/index.html>`__

`SymPy documentation <http://docs.sympy.org/>`__

`SymPy physics vector
documentation <http://docs.sympy.org/latest/modules/physics/vector/index.html>`__

`SymPy mechanics
documentation <http://docs.sympy.org/latest/modules/physics/mechanics/index.html>`__

`PyDy documentation <http://www.pydy.org/documentation.html>`__
