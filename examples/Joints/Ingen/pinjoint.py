from sympy import Symbol
from pydy.bodies import Body
from pydy.joints import PinJoint
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

gravity = Symbol('gravity')
child.add_force((0, child.get_mass() * gravity, 0), child.get_masscenter())

pin_joint = PinJoint('pin_joint', parent, child, parent_point_pos=(1, 1, 1),
                     child_point_pos=(1, 1, 1))

JM = JointsMethod([pin_joint], parent)
print JM.get_kanes()
