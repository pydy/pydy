Once you have SymPy installed on your machine it is time to explore the how it
operates. The first step is to open up your python interpreter. This varies
depending on your work style, but the simplest way is to type python in a
command prompt::

   $ python

Once you are at the python command line the first step is to import basic
functionality from SymPy and the Mechanics module, otherwise you will only have
basic python commands available to work with. We will use the import * method
to bring in all functions from the two modules::

   >>> from sympy import *
   >>> from sympy.physics.mechanics import *

You can now see what functions and variables that are available to you with::

   >>> dir()

This is a long list of available functions, as both packages have many. Read
about the python import statement to learn about better ways to import only
what you need. One good explanation is
<http://effbot.org/zone/import-confusion.htm>.

To get started working with vectors we will need to create a reference frame,
as all vectors need to be defined with respect to a reference frame. If you
know the name of the command that you want to use simply use the builtin help
function to bring up the documentation for the function. In our case we need to
use the ReferenceFrame class::

   >>> help(ReferenceFrame)

Press `q` to return to the command line. Now create an inertial reference frame
called N for Newtonian as was described in the help::

   >>> N = ReferenceFrame('N')

Keep in mind that N is the variable name of which the reference frame named 'N'
is stored. It is important to note that `N` is an object and it has properties
and functions associated with it. To see a list of them type::

   >>> dir(N)

Notice that three of the properties are `x`, `y`, and `z`. These are the
orthonormal unit vectors associated with the reference frame and are the
building blocks for creating vectors. We can create a vector by simply
building a linear combination of the unit vectors::

   >>> v = 1 * N.x + 2 * N.y + 3 * N.z

Now a vector expressed in the N reference frame is stored in the variable `v`.
We can print `v` to the screen by typing::

   >>> print(v)
   N.x + 2*N.y + 3*N.z

The vector `v` can be manipulated as expected. You can multiply and divide them
by scalars::

   >>> 2 * v
   2*N.x + 4*N.y + 6*N.z
   >>> v / 3.0
   0.333333333333333*N.x + 0.666666666666667*N.y + N.z

Note that three is expressed as `3.0` instead of `3`. The python language does
integer division by default. There are ways around this, but for now simply
remember to always declare numbers as floats (i.e. include a decimal).

You can add and subtract vectors::

   >>> v + v
   2*N.x + 4*N.y + 6*N.z
   >>> w = 5 * N.x + 7 * N.y
   >>> v - w
   - 4*N.x - 5*N.y + 3*N.z

Vectors also have some useful properties::

   >>> dir(v)

You can find the magnitude of a vector by typing::

   >>> v.magnitude()
   sqrt(14)

You can compute a unit vector in the direction of `v`::

   >>> v.normalize()
   sqrt(14)/14*N.x + sqrt(14)/7*N.y + 3*sqrt(14)/14*N.z

You can find the measure numbers and the reference frame the vector was defined in
with::

   >>> v.args
   [([1]
   [2]
   [3], N)]

Dot and cross products of vectors can also be computed::

   >>> dot(v, w)
   19
   >>> cross(v, w)
   - 21*N.x + 15*N.y - 3*N.z

We've only used numbers as our measure numbers so far, but it is just as easy
to use symbols. We will introduce six symbols for our measure numbers with the
SymPy `symbols` [`help(symbols) for the documentation`] function::

   >>> a1, a2, a3 = symbols('a1 a2 a3')
   >>> b1, b2, b3 = symbols('b1 b2 b3')

And create two new vectors that are completely symbolic::

   >>> x = a1 * N.x + a2 * N.y + a3 * N.z
   >>> y = b1 * N.x + b2 * N.y + b3 * N.z
   >>> dot(x, y)
   a1*b1 + a2*b2 + a3*b3
   >>> cross(x, y)
   (a2*b3 - a3*b2)*N.x + (-a1*b3 + a3*b1)*N.y + (a1*b2 - a2*b1)*N.z

Numbers and symbols work together seamlessly::

   >>> dot(v, x)
   a1 + 2*a2 + 3*a3

You can also differentiate a vector with respect to a variable in a
reference frame::

   >>> x.diff(a1, N)
   N.x
   >>> z.diff(a1, A)
   - b3*N.y + b2*N.z

Vectors don't have be defined with respect to just one reference frame. We can
create a new reference frame and orient it with respect to the `N` frame that
has already been created. We will use the `orient` method of the new frame to
do a simple rotation through `alpha` about the `N.x` axis::

   >>> A = ReferenceFrame('A')
   >>> alpha = symbols('alpha')
   >>> A.orient(N, 'Axis', [alpha, N.x])

Now the direction cosine matrix with of `A` with respect to `N` can be
computed::

   >>> A.dcm(N)
   [1,           0,          0]
   [0,  cos(alpha), sin(alpha)]
   [0, -sin(alpha), cos(alpha)]

Now that SymPy knows that `A` and `N` are oriented with respect to each other
we can express the vectors that we originally wrote in the `A` frame::

   >>> v.express(A)
   A.x + (3*sin(alpha) + 2*cos(alpha))*A.y + (-2*sin(alpha) + 3*cos(alpha))*A.z
   >>> z = cross(x, y)
   >>> z.express(A)
   >>> (a2*b3 - a3*b2)*A.x + ((a1*b2 - a2*b1)*sin(alpha) + (-a1*b3 +
   a3*b1)*cos(alpha))*A.y + ((a1*b2 - a2*b1)*cos(alpha) + (a1*b3 - a3*b1)*sin(alpha))*A.z

In dynamics systems at least some of the relative orientation of reference
frames and vectors are time varying. The mechanics module provides a way to
specify quantities as time varying. Let's define two variables `beta` and `d` as
variables which are functions of time::

   >>> beta, d = dynamicsymbols('beta d')

Now we can create a new reference frame that is oriented with respect to the A
frame by `beta` and create a vector in that new frame that is a function of `d`.
This time we will use the `orientnew` method of the `A` frame to create the new
reference frame `B`::

   >>> B = A.orientnew('B', 'Axis', [beta, A.y])
   >>> vec = d * B.z

We can now compute the angular velocity of the reference frame `B` with respect
to other reference frames::

   >>> B.ang_vel_in(N)
   beta'*A.y

This allows us to now differentiate the vector, `vec`, with respect to time and
a reference frame::

   >>> vecdot = vec.dt(N)
   >>> vecdot
   d*beta'*B.x + d'*B.z
   >>> vecdot.express(N)
   (d*cos(beta)*beta' + sin(beta)*d')*N.x + (d*sin(A)*sin(beta)*beta' -
   sin(A)*cos(beta)*d')*N.y + (-d*sin(beta)*cos(A)*beta' +
   cos(A)*cos(beta)*d')*N.z

The `dynamicsymbols` function also allows you to store the derivatives of time
varying variables. For example, we can define `omega` as the first time
derivative of `beta` which allows you to explicitly interact with the
derivatives::

   >>> theta = dynamicsymbols('theta')
   >>> omega = dynamicsymbols('theta', 1)

At this point we now have all the tools needed to setup the kinematics for a
dynamic system. Let's start with a simple system where a point can move back
and forth on a spinning disc. First create an inertial reference frame::

   >>> N = ReferenceFrame('N')

Now create a reference frame for the disc::

   >>> D = ReferenceFrame('D')

The disc rotates with respect to `N` about the `N.x` axis through `theta`::

   >>> theta = dynamicsymbols('theta')
   >>> D.orient(N, 'Axis', [theta, N.x])

Define one point at the origin of rotation which is fixed in `N`::

   >>> no = Point('no')
   >>> no.set_vel(N, 0)

The second point, `p`, can move through `x` along the `D.y` axis::

   >>> p = Point('p')
   >>> r = x * D.y
   >>> p.set_pos(no, r)
   >>> p.set_vel(D, r.dt())

The velocity of the point in the `N` frame can now be computed::

   >>> p.vel(N)
   x'*D.y + x*theta'*D.z
   >>> p.vel(N).express(N)
   (-x*sin(theta)*theta' + cos(theta)*x')*N.y + (x*cos(theta)*theta' + sin(theta)*x')*N.z

The acceleration of the point can also be computed, but for this we will make
use of the fact that both `do` and `no` are in the reference frame and use the
`a2pt_theory` method::

   >>> p.a2pt_theory(no, N, D)
   - x*theta'**2*D.y + x*theta''*D.z

Now let's imagine that the point, `p`, is constrained by a spring with
stiffness, `k`, to point, `no`. This means that a force 
