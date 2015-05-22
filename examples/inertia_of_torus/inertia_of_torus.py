"""Derivation of central inertia of a uniform density torus.

See Kane & Levinson, 1985, Chapter 3, Section 3 for further background.


"""
from sympy import integrate, pi, Matrix, symbols
from sympy.physics.mechanics import dot, cross, inertia, ReferenceFrame

phi, theta, s, R, r, m = symbols('phi theta s R r m')

# Volume and density of a torus
V = 2 * R * (pi * r)**2
rho = m / V

A = ReferenceFrame('A')                    # Torus fixed, x-y in symmetry plane
B = A.orientnew('B', 'Axis', [phi, A.z])   # Intermediate frame
C = B.orientnew('C', 'Axis', [-theta, B.y])# Intermediate frame

# Position vector from torus center to arbitrary point of torus
# R : torus major radius
# s : distance >= 0 from center of torus cross section to point in torus
p = R * B.x + s * C.x

# Determinant of the Jacobian of the mapping from a, b, c to x, y, z
# See Wikipedia for a lucid explanation of why we must comput this Jacobian:
# http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant#Further_examples
J = Matrix([dot(p, uv) for uv in A]).transpose().jacobian([phi, theta, s])
dv = J.det().trigsimp()      # Need to ensure this is positive
print("dx*dy*dz = {0}*dphi*dtheta*ds".format(dv))

# We want to compute the inertia scalars of the torus relative to it's mass
# center, for the following six unit vector pairs
unit_vector_pairs = [(A.x, A.x), (A.y, A.y), (A.z, A.z),
                     (A.x, A.y), (A.y, A.z), (A.x, A.z)]

# Calculate the six unique inertia scalars using equation 3.3.9 of Kane &
# Levinson, 1985.
inertia_scalars = []
for n_a, n_b in unit_vector_pairs:
    # Integrand of Equation 3.3.9
    integrand = rho * dot(cross(p, n_a), cross(p, n_b)) * dv

    # Compute the integral by integrating over the whole volume of the tours
    I_ab = integrate(integrate(integrate(integrand,
                 (phi, 0, 2*pi)), (theta, 0, 2*pi)), (s, 0, r))

    inertia_scalars.append(I_ab)

# Create an inertia dyad from the list of inertia scalars
I_A_O = inertia(A, *inertia_scalars)
print("Inertia of torus about mass center = {0}".format(I_A_O))
