from sympy import Symbol
from pydy.bodies import Body
from pydy.joints import SlidingJoint
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

gravity = Symbol('gravity')
child.add_force((child.get_mass() * gravity, 0, 0), child.get_masscenter())

sliding_joint = SlidingJoint('sliding_joint', parent, child,
                             parent_point_pos=(1, 1, 1),
                             child_point_pos=(1, 1, 1))

spring_constant = Symbol('k')
child.add_force((spring_constant * sliding_joint.get_coordinates()[0], 0, 0),
                child.get_masscenter())

JM = JointsMethod([sliding_joint], parent)
print JM.get_kanes()
