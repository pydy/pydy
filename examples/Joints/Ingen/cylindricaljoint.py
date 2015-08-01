from sympy import Symbol
from pydy.bodies import Body
from pydy.joints import CylindricalJoint
from pydy.joints_method import JointsMethod

parent = Body('parent')
child = Body('child')

cylindrical_joint = CylindricalJoint('cylindrical_joint')