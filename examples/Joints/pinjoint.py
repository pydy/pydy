from sympy import symbols, Symbol, acos
from sympy.physics.mechanics import *
from sympy.physics.vector import dot

from numpy import linspace
from pydy.system import System
from numpy import sqrt

from pydy.viz.shapes import Cylinder, Sphere
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

# Body
# -----------------------------------------------------------------------------
# Ground/ Root Body
root_frame = ReferenceFrame('root_frame')
origin = Point('origin')
origin.set_vel(root_frame, 0)

# Parent Body
parent_frame = ReferenceFrame('parent_frame')
parent_frame.set_ang_vel(root_frame, 0)
parent_masscenter = Point('parent_masscenter')
parent_mass = Symbol('parent_mass')
parent = Particle('ParP', parent_masscenter, parent_mass)


# Child Body
child_frame = ReferenceFrame('child_frame')
child_masscenter = Point('child_masscenter')
child_mass = Symbol('child_mass')
child = Particle('ParR', child_masscenter, child_mass)

# Connecting ground and parent (without a Joint)
# -----------------------------------------------------------------------------
parent_frame.orient(root_frame, 'Axis', [0, root_frame.z])
# Here, parent is at origin, but user will be able to add a pinjoint between
# origin i.e. Ground and parent too.
parent_masscenter.set_pos(origin, 0)
parent_masscenter.v2pt_theory(origin, root_frame, parent_frame)

# PinJoint between parent and child.
# -----------------------------------------------------------------------------
theta = dynamicsymbols('theta')
thetad = dynamicsymbols('theta', 1)
omega = dynamicsymbols('omega')
omegad = dynamicsymbols('omega', 1)
g = symbols('g')

child_frame.orient(root_frame, 'Axis', [theta, root_frame.z])
child_frame.set_ang_vel(root_frame, omega * root_frame.z)

parent_joint_point = parent_masscenter.locatenew(
    'parent_joint',
    parent_frame.x + parent_frame.y + parent_frame.z)

child_joint_point = child_masscenter.locatenew(
    'child_joint',
    child_frame.x + child_frame.y + child_frame.z)

child_joint_point.set_pos(parent_joint_point, 0)
child_masscenter.v2pt_theory(parent_masscenter, root_frame, child_frame)

# JointsMethod generating equations of motion
# -----------------------------------------------------------------------------

kd = [thetad - omega]
FL = [(child_masscenter, child_mass * g * child_frame.y)]
BL = [parent, child]

KM = KanesMethod(root_frame, q_ind=[theta], u_ind=[omega], kd_eqs=kd)
print KM.kanes_equations(FL, BL)

# Numerical part
# -----------------------------------------------------------------------------

constants = {child_mass: 10.0, g: 9.81}

initial_conditions = {theta: 1.0, omega: 0.0}

sys = System(KM, constants=constants, initial_conditions=initial_conditions)

frames_per_sec = 60
final_time = 5.0

times = linspace(0.0, final_time, final_time * frames_per_sec)
sys.times = times
x = sys.integrate()

# Visualization part
# -----------------------------------------------------------------------------
link1 = Cylinder(name='link1', radius=0.1,
                 length=sqrt(3), color='red')
link2 = Cylinder(name='link2', radius=0.1,
                 length=sqrt(3), color='red')
joint_link = Cylinder(name='joint_link', radius=0.05, length=0.5, color='blue')
sphere = Sphere(name='sphere', radius=0.2, color='green')
joint_link_sphere = Sphere(name='sphere', radius=0.1, color='blue')

# parent
sphereP_viz_frame = VisualizationFrame('sphereP', root_frame, parent_masscenter,
                                       sphere)

# parent to joint link
parent_joint_vec = parent_joint_point.pos_from(parent_masscenter)
parent_angle = acos(dot(parent_frame.y, parent_joint_vec)/parent_joint_vec.magnitude())
linkR_frame = parent_frame.orientnew('frameR', 'Axis', [-parent_angle, cross(parent_joint_vec, root_frame.y)])
linkR_link_half = parent_masscenter.locatenew('originP', parent_joint_vec/2)
link_parent_joint = VisualizationFrame('linkP_joint', linkR_frame, linkR_link_half, link1)

# joint to child link
child_joint_vec = child_masscenter.pos_from(child_joint_point)
child_angle = acos(dot(child_frame.y, child_joint_vec)/child_joint_vec.magnitude())
tmp = cross(child_joint_vec, child_frame.y)
linkJ_frame = child_frame.orientnew('frameJ', 'Axis', [-child_angle, tmp])
child_joint_vec = child_masscenter.pos_from(child_joint_point)
linkJ_link_half = parent_joint_point.locatenew('originJ', child_joint_vec/2)
link_joint_child = VisualizationFrame('linkJ_joint', linkJ_frame, linkJ_link_half, link1)

# joint link
joint_link_vec = cross(parent_joint_vec, child_joint_vec)
joint_link_angle = acos(dot(parent_frame.y, joint_link_vec)/joint_link_vec.magnitude())
link_joint_frame = child_frame.orientnew('frame_joint', 'Axis', [joint_link_angle, child_frame.y])
link_joint = VisualizationFrame('joint_link', root_frame, parent_joint_point,
                                joint_link_sphere)

# child
sphereR_viz_frame = VisualizationFrame('sphereR', root_frame, child_masscenter,
                                       sphere)

world_frame = root_frame.orientnew('world', 'Axis', [0, root_frame.z])
scene = Scene(world_frame, origin, link_parent_joint, link_joint, link_joint_child,
              sphereP_viz_frame, sphereR_viz_frame)

scene.generate_visualization_json_system(sys)
scene.display()
