from sympy import symbols, sin, cos, tan
from sympy.physics.mechanics import *

# -- Declare symbols.

# Declare degrees of freedom, and their time derivatives.
# TODO get rid of those which are unnecessary.
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
frameA = N.orientnew('frameA', 'Axis', [qA, N.x])
frameB = frameA.orientnew('frameB', 'Axis', [qB, frameA.z])

# TODO why do we need to do this?
frameA.set_ang_vel(N, uA * N.x)
frameA.set_ang_acc(N, frameA.ang_vel_in(N).dt(N)) # TODO

frameB.set_ang_vel(frameA, uB * frameA.z)
frameB.set_ang_acc(frameA, frameB.ang_vel_in(frameA).dt(frameA))

# Origin.
NO = Point('NO')
NO.set_vel(N, 0)

# Center of mass of upper rod.
Acm = NO.locatenew('Acm', -LA * frameA.z)
Acm.v2pt_theory(NO, N, frameA) # Don't need this in MotionGenesis.
Bcm = NO.locatenew('Bcm', - LB * frameA.z)
Bcm.v2pt_theory(NO, N, frameA) # Don't need this in MotionGenesis.

# Inertia dyadic.
# TODO inertia dyadics are about a specific point, right? are we requiring the user to define it about the COM always or something? idk if this info should be separate. we have inertia info floating about in 2 different locations.
IA = inertia(frameA, IAx, 0, 0)
IB = inertia(frameB, IBx, IBy, IBz)

# Create rigid bodies.
BodyA = RigidBody('BodyA', Acm, frameA, mA, (IA, Acm))
BodyB = RigidBody('BodyB', Bcm, frameB, mB, (IB, Bcm))
BodyList = [BodyA, BodyB]

# Forces.
# Would be nice to have a method that applies gravity force to all objects.
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
