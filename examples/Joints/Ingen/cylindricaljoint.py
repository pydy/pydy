from pydy.bodies import Body
from pydy.joints import CylindricalJoint
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

cylindrical_joint = CylindricalJoint('cylindrical_joint',parent, child,
                                     child_point_pos=(1, 1, 1))

JM = JointsMethod([cylindrical_joint], parent)
print JM.get_kanes()
