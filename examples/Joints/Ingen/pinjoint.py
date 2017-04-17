from numpy import linspace

from sympy import Symbol
from pydy.bodies import Body
from pydy.joints import PinJoint
from pydy.system import System
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

gravity = Symbol('gravity')
child.add_force((0, child.get_mass() * gravity, 0), child.get_masscenter())

pin_joint = PinJoint('pin_joint', parent, child, parent_point_pos=(1, 1, 1),
                     child_point_pos=(1, 1, 1))

JM = JointsMethod([pin_joint], parent)

# Sample implementation since numeric part is not yet implemented.
# -----------------------------------------------------------------------------
constants = {child.get_mass(): 10.0}

theta = pin_joint.get_coordinates()[0]
omega = pin_joint.get_speeds()[0]
initial_conditions = {theta: 1.0, omega: 0.0}

sys = System(JM, constants=constants, initial_conditions=initial_conditions)

frames_per_sec = 60
final_time = 5.0

times = linspace(0.0, final_time, final_time * frames_per_sec)
sys.times = times
x = sys.integrate()
