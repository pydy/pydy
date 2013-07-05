from sympy.physics.mechanics import *
#from pydy-viz import Shape, Cube, 

def test_shape():
    point_list = [(1, 0, 0), (2, 3, 1), (-1, 5, -1), (-4, 2, 5)]
    shape = Shape('shape')
    
    #checking setter methods ..
    shape.set_point_list(point_list)
    shape.set_color('RED')  #abstract method (abc package)
    shape.set_center([0, 0, 0]) #abstract method
    
    #checking getter methods ...
    assert shape.get_name() == 'shape' #abstract method
    
    #How to compare lists?
    assert shape.get_point_list() ==  [(1, 0, 0), (2, 3, 1), \
                                      (-1, 5, -1), (-4, 2, 5)]
    
    assert shape.get_color() == 'RED' #abstract method
    #Again how to compare tuples
    assert shape.get_center() == [0, 0, 0]  #abstractmethod


def test_cube():
    cube = Cube('cube')
    cube.set_length(10)
    cube.set_color('BLUE')
    cube.set_center([0, 0, 0])
    
    assert cube.get_name() == 'Cube'
    assert cube.get_length() == 10
    assert cube.get_point_list() ==  [(5, 5, 5), (5, 5, -5), \
                                     (5, -5, 5), (-5, -5, -5), \
                                       (-5 ,5 , 5), (-5, 5, -5), \
                                          (-5, -5, 5), (-5, -5, -5)]
    assert cube.color() == 'BLUE'
    assert cube.get_center() == [0, 0, 0]
    
    assert isinstance(cube, Shape) == True

def test_cylinder():
    cylinder = Cylinder('cylinder')
    cylinder.set_height(10)
    cylinder.set_radius(5)
    cylinder.set_color('BLACK')
    cylinder.set_center([0, 0, 0])
    
    
    assert cylinder.get_name() == 'cylinder'
    assert cylinder.get_height() == 10
    assert cylinder.get_radius() == 5
    assert cylinder.color() == 'BLACK'
    assert cylinder.get_center() == (0, 0, 0)
    
    assert isinstance(cylinder,Shape) == True
    
def test_cone():
    cone = Cone('cone')
    cone.set_height(10)
    cone.set_radius(5)
    cone.set_color('GREEN')
    cone.set_center([0, 0, 0])	 #base center, or geometric center?
    
    assert cone.get_name() == 'cone'
    assert cone.get_height() == 10
    assert cone.get_radius() == 5
    assert cone.get_color() == 'GREEN'                                            
    assert cone.get_center() == [0, 0, 0]
    
    assert isinstance(cone,Shape) == True
                            
	                                                       
 

