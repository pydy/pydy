import sympy as sm
from sympy.physics.vector import *
from sympy.physics.mechanics import *
from numpy import linspace
from pydy.system import System
from numpy import pi
from pydy.viz.shapes import Cylinder, Sphere
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

# Body----------------------------------------------
parent_masscenter = Point('parent_masscenter')
parent_frame = ReferenceFrame('parent_frame')
parent_mass = sm.Symbol('parent_mass')

parent_masscenter.set_vel(parent_frame, 0)
parent = Particle('parent', parent_masscenter, parent_mass)

child_masscenter = Point('child_masscenter')
child_frame = ReferenceFrame('child_frame')
child_mass = sm.Symbol('child_mass')

child_masscenter.set_vel(child_frame, 0)
child = Particle('child', child_masscenter, child_mass)

print "Finished creating Bodies"
# Joint----------------------------------------------

# join parent and child frames
child_frame.orient(parent_frame, 'Axis', [0, parent_frame.x])

coordinate = dynamicsymbols('theta')
coordinated = dynamicsymbols('theta', 1)
omega = dynamicsymbols('omega')

a, b, c, d, e, f, g = sm.symbols('a b c d e f g')

child_joint_loc = (a, b, c)
parent_joint_loc = (d, e, f)

parent_joint_point = parent.get_point().locatenew(
    'parent_joint',
    a * parent_frame.x + b * parent_frame.y + c * parent_frame.z)

child_joint_point = child.get_point().locatenew(
    'child_joint',
    e * child_frame.x + d * child_frame.y + f * child_frame.z)

child_joint_point.set_pos(parent_joint_point, 0)


child_frame.orient(parent_frame, 'Axis', [coordinate, parent_frame.x])
child_frame.set_ang_vel(parent_frame, omega * parent_frame.x)

child.get_point().v2pt_theory(parent.get_point(), parent_frame, child_frame)

print "Finished creating Joints"
# JointsMethod----------------------------------------------------
BL = [parent, child]
FL = [(parent.get_point(), parent.get_mass() * g * parent_frame.x), (child.get_point(), child.get_mass() * g * parent_frame.x)]
kd = [coordinated - omega]

KM = KanesMethod(parent_frame, q_ind=[coordinate], u_ind=[omega], kd_eqs=kd)
(fr, frstar) = KM.kanes_equations(FL, BL)
print "###########", fr
print "###########", frstar

print "Finished JointsMethod"
# Numeric part---------------------------------------------------
constants = {parent_mass: 1.0, child_mass: 1.0, a: 1.0, b: 1.0, c: 1.0, d: 1.0, e: 1.0, f: 1.0, g: 9.8}
#constants = {g: 9.8}
initial_conditions = {coordinate: 0.0, omega: 0.0}
sys = System(KM, constants=constants, initial_conditions=initial_conditions)

frames_per_sec = 60
final_time = 5

times = linspace(0.0, final_time, final_time * frames_per_sec)
sys.times = times
x = sys.integrate()

print "Numeric part"
# Visualization part-------------------------------------------
parent_sphere = Sphere(name='sphere', radius=1.0)
child_sphere = Sphere(name='sphere', radius=1.0)

sphereP_viz_frame = VisualizationFrame('sphereP', parent_frame, parent_masscenter, parent_sphere)
sphereR_viz_frame = VisualizationFrame('sphereR', child_frame, child_masscenter, child_sphere)

scene = Scene(parent_frame, parent_masscenter, sphereP_viz_frame, sphereR_viz_frame)
scene.generate_visualization_json_system(sys)

scene.display()

