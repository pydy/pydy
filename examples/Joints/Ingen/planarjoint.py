from sympy import Symbol
from pydy.bodies import Body
from pydy.joints import PlanarJoint
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

planar_joint = PlanarJoint('planar_joint', parent, child, (0, 0, 0), (0, 1, 0),
                           (0, 1, 0), (0, 1, 0))

k = Symbol('k')
disx = planar_joint.get_coordinates()[1]
disy = planar_joint.get_coordinates()[2]

child.add_force((k * disx, 0, 0), child.get_masscenter())
child.add_force((0, k * disy, 0), child.get_masscenter())

JM = JointsMethod([planar_joint], parent)
print JM.get_kanes()
