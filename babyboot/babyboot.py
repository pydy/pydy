from sympy import symbols, sin, cos, tan
from sympy.physics.mechanics import *

# -- Declare symbols.

# Declare degrees of freedom, and their time derivatives.
qA, qB, uA, uB = dynamicsymbols('qA qB uA uB')
qAd, qBd, uAd, uBd = dynamicsymbols('qA qB uA uB', 1)

# Declare the constants in the problem.
LA, LB, mA, mB, IAx, IBx, IBy, IBz, g = symbols(
        'LA LB mA mB IAx IBx IBy IBz g')

# TODO no clue.
mechanics_printing()

# -- Set up geometry.
# Fixed frame at base of upper rod.
N = ReferenceFrame('N')
# Frame tracking the upper rod.
Aframe = N.orientnew('Aframe', 'Axis', [qA, N.x])
Bframe = Aframe.orientnew('Bframe', 'Axis', [qB, Aframe.z])

# TODO why do we need to do this?
Aframe.set_ang_vel(N, qA * N.x)
Aframe.set_ang_acc(N, Aframe.ang_vel_in(N).dt(Aframe)) # TODO

# Origin.
O = Point('O')
O.set_vel(N, 0)

# Center of mass of upper rod.
Acm = O.locatenew('Acm', - LA * Aframe.z)
Acm.set_vel(N, LA * uA * Aframe.y)  # Don't need this in MotionGenesis.
Bcm = O.locatenew('Bcm', - LB * Aframe.z)
Bcm.set_vel(N, LB * uA * Aframe.y) # Don't need this in MotionGenesis.

# Inertia dyadic.
IA = inertia(Aframe, IAx, 0, 0)
IB = inertia(Bframe, IBx, IBy, IBz)

# Create rigid bodies.
BodyA = RigidBody('BodyA', Acm, Aframe, mA, (IA, Acm))
BodyB = RigidBody('BodyB', Bcm, Bframe, mB, (IB, Bcm))
BodyList = [BodyA, BodyB]

# Forces.
ForceList = [(Acm, - mA * g * N.z), (Bcm, - mB * g * N.z)]

# Kinematic differential equations. TODO necessary?
kd = [qAd -uA, qBd - uB] # TODO constrain upper rod to O?

KM = KanesMethod(N, q_ind=[qA, qB], u_ind=[uA, uB], kd_eqs=kd)
(fr, frstar) = KM.kanes_equations(ForceList, BodyList)

# Get equations of motion.
MM = KM.mass_matrix
forcing = KM.forcing
rhs = MM.inv() * forcing
kdd = KM.kindiffdict()
rhs = rhs.subs(kdd)
rhs.simplify()
mprint(rhs)
