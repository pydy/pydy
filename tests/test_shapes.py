from sympy.physics.mechanics import *
#from pydy-viz import Shape,Cube

def test_shape():
    point_list = [(1, 0, 0), (2, 3, 1), (-1, 5, -1), (-4, 2, 5)]
    shape = Shape('shape', pt_list=point_list, color="RED", \
                                            origin=(0, 0, 0))
    assert shape.get_color() == 'RED'
    #How to compare lists?
    assert shape.get_pt_list() ==  [(1, 0, 0), (2, 3, 1), \
                                      (-1, 5, -1), (-4, 2, 5)]
    assert shape.get_name() == 'shape'


def test_cube():
    cube = Cube('cube',length=10,color="BLUE",origin=(0,0,0))
    assert cube.get_color() == 'BLUE'
    assert cube.get_pt_list() ==  [(5, 5, 5), (5, 5, -5), \
                                     (5, -5, 5), (-5, -5, -5), \
                                       (-5 ,5 , 5), (-5, 5, -5), \
                                          (-5, -5, 5), (-5, -5, -5)]
    assert cube.get_name() == 'Cube'
    assert isinstance(cube,Shape) == True

#def test_cylinder():
 

