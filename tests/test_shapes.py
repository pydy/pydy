from sympy.physics.mechanics import *
from numpy.testing import assert_allclose
#from pydy-viz import Shape, Cube, 

def test_shape():
    point_list = [[1, 0, 0], [2, 3, 1], [-1, 5, -1], [-4, 2, 5]]
    shape = Shape('shape')
    
    shape.point_list = point_list
    shape.color = 'RED'  #abstract method (abc package), a @property  
    shape.center = [0, 0, 0] #abstract method
    
    assert shape.name == 'shape' #abstract method, a @property
    assert_allclose(shape.point_list, point_list)
    assert shape.color == 'RED' #abstract method
    assert_allclose(shape.center, [0, 0, 0])  


def test_cube():
    cube = Cube('cube')
    cube.length = 10
    cube.color = 'BLUE'
    cube.center = [0, 0, 0]
    
    assert cube.name == 'Cube'
    assert cube.length == 10
    assert cube.color == 'BLUE'
    point_list = [[5, 5, 5], [5, 5, -5], \
                    [5, -5, 5], [-5, -5, -5], \
                        [-5, 5, 5], [-5, 5, -5], \
                           [-5, -5, 5], [-5, -5, -5]]

    assert_allclose(cube.point_list,point_list)
    assert_allclose(cube.center, [0, 0, 0])
    
    assert isinstance(cube, Shape) == True

def test_cylinder():
    cylinder = Cylinder('cylinder')
    cylinder.height = 10
    cylinder.radius = 5
    cylinder.color = 'BLACK'
    cylinder.center = [0, 0, 0]
    
    
    assert cylinder.name == 'cylinder'
    assert cylinder.height == 10
    assert cylinder.radius == 5
    assert cylinder.color == 'BLACK'
    assert_allclose(cylinder.center, [0, 0, 0])
    
    assert isinstance(cylinder, Shape) == True
    
def test_cone():
    cone = Cone('cone')
    cone.height = 10
    cone.radius = 5
    cone.color = 'GREEN'
    cone.center = [0, 0, 0]	 #base center, or geometric center?
    
    assert cone.name == 'cone'
    assert cone.height == 10
    assert cone.radius == 5
    assert cone.color == 'GREEN'                                            
    assert_allclose(cone.center, [0, 0, 0])
    
    assert isinstance(cone, Shape) == True
                            
def test_sphere():
    sphere = Sphere('sphere')
    sphere.radius = 10
    sphere.color = 'GREY'
    sphere.center = [0, 0, 0]	                                                       
    
    assert sphere.name == 'sphere'
    assert sphere.radius == 10
    assert sphere.color == 'GREY'
    assert_allclose(sphere.center, [0, 0, 0])
  
    assert isinstance(sphere, Shape) == True
 


def test_circle():
    circle = Circle('circle')
    circle.radius = 10
    circle.color = 'CYAN'
    circle.center = [0, 0, 0]

    assert circle.name == 'circle'
    assert circle.radius == 10
    assert circle.color == 'CYAN'
    assert_allclose(circle.center, [0, 0, 0])

    assert isinstance(circle, Shape) == True

def test_mesh_shape():
    mesh_shape = MeshShape('mesh_shape')
    point_list = [[2, 3, 1], [4, 6, 2], \
                         [5, 3, 1], [5, 3, 6], \
                         [2, 8, 4], [7, 4, 1]]
    mesh_shape.points = point_list
    mesh_shape.color = 'PINK'
    mesh_shape.origin = [0, 0, 0]
    
    assert mesh_shape.name == 'mesh_shape'
    assert_allclose(mesh_shape.points,point_list)
    assert mesh_shape.color == 'PINK'
    assert_allclose(mesh_shape.origin, point_list)
    
    assert isinstance(mesh_shape, Shape)                      

def test_tetrahedron():
    pass

def test_octahedron():
    pass

def test_icosahedron():
    pass


def test_plane():
    plane = Plane('plane')
    plane.length = 10
    plane.width = 20
    plane.color = 'BLACK'
    plane.center = [0, 0, 0]
    
    
    assert plane.name == 'plane'
    assert plane.length == 10
    assert plane.width == 20
    assert plane.color == 'BLACK'
    assert_allclose(plane.center, [0, 0, 0])
    
    assert isinstance(plane, Shape) == True
 
def test_torus():
    pass

def test_tube():
    pass

def test_torus_knot():
    pass

 
