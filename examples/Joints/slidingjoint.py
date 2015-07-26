from sympy import Symbol, acos, sin
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

# PinJoint between parent and child.
# -----------------------------------------------------------------------------
parent_masscenter = Point('parent_masscenter')
parent_mass = Symbol('parent_mass')
parent = Particle('ParP', parent_masscenter, parent_mass)


# Child Body
child_frame = ReferenceFrame('child_frame')
child_masscenter = Point('child_masscenter')
child_mass = Symbol('child_mass')
child_masscenter.set_vel(child_frame, 0)
child = Particle('ParR', child_masscenter, child_mass)

# Connecting ground and parent (without a Joint)
# -----------------------------------------------------------------------------
parent_frame.orient(root_frame, 'Axis', [0, root_frame.z])
# Here, parent is at origin, but user will be able to add a pinjoint between
# origin i.e. Ground and parent too.
parent_masscenter.set_pos(origin, 0)
parent_masscenter.v2pt_theory(origin, root_frame, parent_frame)


# SlidingJoint between parent and child.
# -----------------------------------------------------------------------------
dis = dynamicsymbols('dis')
disd = dynamicsymbols('dis', 1)
vel = dynamicsymbols('vel')
veld = dynamicsymbols('vel', 1)
g = Symbol('g')

child_frame.orient(parent_frame, 'Axis', [0, root_frame.z])

parent_joint_point = parent_masscenter.locatenew(
    'parent_joint',
    parent_frame.x + parent_frame.y + parent_frame.z)

child_joint_point = child_masscenter.locatenew(
    'child_joint',
    - child_frame.x)

parent_joint_point.set_vel(parent_frame, 0)
child_joint_point.set_vel(child_frame, 0)

direction = parent_frame.x

child_joint_point.set_pos(parent_joint_point, dis * direction)
child_joint_point.set_vel(parent_frame, vel * direction)
child_masscenter.set_vel(root_frame, vel*direction)

# JointsMethod generation equations of motion
# -----------------------------------------------------------------------------

kd = [disd - vel]
FL = [(child_masscenter, child_mass * g * root_frame.x),
      (child_masscenter, -0.001 * dis * root_frame.x)]
BL = [parent, child]

KM = KanesMethod(root_frame, q_ind=[dis], u_ind=[vel], kd_eqs=kd)
(fr, frstar) = KM.kanes_equations(FL, BL)
print fr
print frstar

# Numerical part
#  ---------------------------------------------------------------------------
constants = {child_mass: 0.1, g: 0.001}

initial_conditions = {vel: 0.01, dis: 0}

sys = System(KM, constants=constants, initial_conditions=initial_conditions)

frames_per_sec = 60
final_time = 3.0

times = linspace(0.0, final_time * frames_per_sec)
sys.times = times
x = sys.integrate()

# Visualization part
# -----------------------------------------------------------------------------
body_sphere = Sphere(name='body_sphere', radius=0.2, color='green')
joint_sphere = Sphere(name='body_sphere', radius=0.15, color='yellow')
joint_sphere2 = Sphere(name='body_sphere', radius=0.15, color='grey')
link = Cylinder(name='link', radius=0.1, length=sqrt(3), color='red')
link2 = Cylinder(name='link', radius=0.1, length=sqrt(3)/2, color='red')
joint_link = Cylinder(name='joint_link', radius=0.05, length=0.5, color='blue')

parent_viz_frame = VisualizationFrame('sphereP', root_frame, parent_masscenter,
                                      body_sphere)

child_viz_frame = VisualizationFrame('sphereC', root_frame, child_masscenter,
                                     body_sphere)

childjoint_viz_frame = VisualizationFrame('sphereCJ', root_frame, child_joint_point,
                                          joint_sphere)
parentjoint_viz_frame = VisualizationFrame('spherePJ', root_frame, parent_joint_point,
                                           joint_sphere2)


# parent to joint link
parent_joint_vec = parent_joint_point.pos_from(parent_masscenter)
parent_angle = acos(dot(parent_frame.y, parent_joint_vec)/parent_joint_vec.magnitude())
linkR_frame = parent_frame.orientnew('frameR', 'Axis',
                                     [-parent_angle, cross(parent_joint_vec, root_frame.y)])
linkR_link_half = parent_masscenter.locatenew('originP', parent_joint_vec/2)
link_parent_joint = VisualizationFrame('linkP_joint', linkR_frame, linkR_link_half, link)

# child to joint link
child_joint_vec = child_masscenter.pos_from(child_joint_point)
child_angle = acos(dot(child_frame.y, child_joint_vec)/child_joint_vec.magnitude())
tmp = cross(child_joint_vec, child_frame.y)
linkJ_frame = child_frame.orientnew('frameJ', 'Axis', [-child_angle, tmp])
linkJ_link_half = child_joint_point.locatenew('originJ', child_joint_vec/2)
link_joint_child = VisualizationFrame('linkJ_joint', linkJ_frame, linkJ_link_half,
                                      link2)

# joint
joint_link_vec = cross(parent_joint_vec, child_joint_vec)
joint_link_angle = acos(dot(parent_frame.y, joint_link_vec)/joint_link_vec.magnitude())
link_joint_frame = child_frame.orientnew('frame_joint', 'Axis',
                                         [joint_link_angle, root_frame.x])
link_joint = VisualizationFrame('joint_link', link_joint_frame, parent_joint_point,
                                joint_link)

world_frame = root_frame.orientnew('world', 'Axis', [0, root_frame.z])
scene = Scene(world_frame, origin, link_parent_joint, link_joint, link_joint_child,
              parent_viz_frame, child_viz_frame, childjoint_viz_frame, parentjoint_viz_frame)

scene.generate_visualization_json_system(sys)
scene.display()
