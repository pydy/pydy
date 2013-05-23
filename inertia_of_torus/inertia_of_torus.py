"""Derivation of central inertia of a uniform density torus."""
from sympy import integrate, pi, Matrix, symbols
from sympy.physics.mechanics import dot, cross, inertia, Point, ReferenceFrame

phi, theta, s, R, r, m = symbols('phi theta s R r m')

# Volume and density of a torus
V = 2 * R * (pi * r)**2
rho = m / V

# A is the frame fixed in the torus with A.x and A.y in the symmetry plane
A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [phi, A.z])
C = B.orientnew('C', 'Axis', [-theta, B.y])

# Position vector from torus center to arbitrary point of torus
p = R * B.x + s * C.x

# Determinant of the Jacobian of the mapping from a, b, c to x, y, z
J = Matrix([dot(p, uv) for uv in A]).transpose().jacobian([phi, theta, s])
dv = J.det().trigsimp()      # Need to ensure this is positive
print('dx*dy*dz = {0}*dphi*dtheta*ds'.format(dv))

# Integrands of definition of inertia scalars
unit_vector_pairs = [(A.x, A.x), (A.y, A.y), (A.z, A.z),
		     (A.x, A.y), (A.y, A.z), (A.x, A.z)]

inertia_scalars = []
for uv1, uv2 in unit_vector_pairs:
    integrand = rho * dot(cross(p, uv1), cross(p, uv2)) * dv
    integral = integrate(integrate(integrate(integrand,
	       (phi, 0, 2*pi)), (theta, 0, 2*pi)), (s, 0, r))
    
    inertia_scalars.append(integral)
 

I_A_O = inertia(A, *inertia_scalars)
print("Inertia of torus about mass center = {0}".format(I_A_O))

