from sympy.physics.mechanics import *
from numpy.testing import assert_allclose
#from pydy-viz import Shape, Cube, 

def test_cube():
    cube = Cube('cube', length=10, \
                         color='BLUE', center=[0,0,0])
    assert cube.name == 'cube'
    assert cube.length == 10
    assert cube.color == 'BLUE'
    assert_allclose(cube.center, [0, 0, 0])
    
    cube.name = 'cube1'
    assert cube.name == 'cube1'

    cube.length = 16
    assert cube.length == 16

    cube.color = 'RED'
    assert cube.color == 'RED'

    cube.center = [6, 5, 3]
    assert_allclose(cube.center, [6, 5, 3])
    
    assert isinstance(cube, Shape)

def test_cylinder():
    cylinder = Cylinder('cylinder', height=10, radius=5, \
                                    color='BLACK', center=[0,0,0])
    assert cylinder.name == 'cylinder'
    assert cylinder.height == 10
    assert cylinder.radius == 5
    assert cylinder.color == 'BLACK'
    assert_allclose(cylinder.center, [0, 0, 0])
    
    
    cylinder.name = 'cylinder1'
    assert cylinder.name == 'cylinder1'

    cylinder.height = 14
    assert cylinder.height == 14

    cylinder.radius = 7
    assert cylinder.radius == 7

    cylinder.color = 'PINK'
    assert cylinder.color == 'PINK'

    cylinder.center = [3, 6, -2]
    assert_allclose(cylinder.center, [3, 6, -2])
    
    assert isinstance(cylinder, Shape)
    
def test_cone():
    cone = Cone('cone',heigth=10, radius=5, \
                        color='GREEN', center=[0,0,0])
    assert cone.name == 'cone'
    assert cone.height == 10
    assert cone.radius == 5
    assert cone.color == 'GREEN'
    assert_allclose(cone.center, [0, 0, 0])
    
    cone.name = 'cone1'
    assert cone.name == 'cone1'

    cone.height = 16
    assert cone.height == 16

    cone.radius = 3
    assert cone.radius == 3

    cone.color = 'BLUE'
    assert cone.color == 'BLUE'                                            
    
    cone.center = [3, 7, 3]
    assert_allclose(cone.center, [3, 7, 3])
    
    assert isinstance(cone, Shape)
                            
def test_sphere():
    sphere = Sphere('sphere', radius=10, color='GREY', \
                                   center=[0,0,0])
    assert sphere.name == 'sphere1'
    assert sphere.radius == 10
    assert sphere.color == 'GREY'
    assert_allclose(sphere.center, [0, 0, 0])	                                                       
    sphere.name = 'sphere1'
    assert sphere.name == 'sphere1'
   
    sphere.radius = 14
    assert sphere.radius == 14

    sphere.color = 'GREEN'
    assert sphere.color == 'GREEN'

    sphere.center = [5, 2, 6]
    assert_allclose(sphere.center, [5, 2, 6])
  
    assert isinstance(sphere, Shape)
 


def test_circle():
    circle = Circle('circle', radius=10, \
                           color='CYAN',center=[0,0,0])

    assert circle.name == 'circle'
    assert circle.radius == 10
    assert circle.color == 'CYAN'
    assert_allclose(circle.center, [0, 0, 0])

    circle.name = 'circle1'
    assert circle.name == 'circle1'

    circle.radius = 12
    assert circle.radius == 12

    circle.color = 'GREEN'
    assert circle.color == 'GREEN'

    circle.center = [0, -3, 5]
    assert_allclose(circle.center, [0, -3, 5])
    
    assert isinstance(circle, Shape)

def test_mesh_shape():
    point_list = [[2, 3, 1], [4, 6, 2], \
                         [5, 3, 1], [5, 3, 6], \
                         [2, 8, 4], [7, 4, 1]]
    
    mesh_shape = MeshShape('mesh_shape', points=point_list, \
                               color='PINK', origin=[0,0,0])                     
    assert mesh_shape.name == 'mesh_shape'
    assert_allclose(mesh_shape.points,point_list)
    assert mesh_shape.color == 'PINK'
    assert_allclose(mesh_shape.origin, point_list)
    
    mesh_shape.name = 'mesh_shape1'
    assert mesh_shape.name == 'mesh_shape1'
    
    new_point_list = [[3, 4, 12], [2, 4, 4], [3, 2, 41], [2, 5, 4]]
    mesh_shape.points = new_point_list
    assert_allclose(mesh_shape.points, new_point_list)
    
    mesh_shape.color = 'BLACK'
    assert mesh_shape.color == 'BLACK'
    
    mesh_shape.origin = [-5, 1, 2]
    assert_allclose(mesh_shape.origin, [-5, 1, 2])
    
    assert isinstance(mesh_shape, Shape)                      

def test_plane():
    plane = Plane('plane', length=10, width=20, \
                       color='BLUE', center=[0,0,0])
    assert plane.name == 'plane'
    assert plane.length == 10
    assert plane.width == 20
    assert plane.color == 'BLUE'
    assert_allclose(plane.center, [0, 0, 0])
    
    plane.name = 'plane1'
    assert plane.name == 'plane1'
    
    plane.length = 30
    assert plane.length == 30
    
    plane.width = 10
    assert plane.width == 10
    
    plane.color = 'BLACK'
    assert plane.color == 'BLACK'
    
    plane.center = [7, 4, 5]
    assert_allclose(plane.center, [7, 4, 5])
    
    assert isinstance(plane, Shape)

def test_tetrahedron():
    pass

def test_octahedron():
    pass

def test_icosahedron():
    pass


 
def test_torus():
    torus = Torus('torus', radius=10, tube_radius=2, \
                        color='BLUE', center=[0,0,0])

    assert torus.name == 'torus' 			
    assert torus.radius == 10
    assert torus.tube_radius == 2
    assert torus.color == 'BLUE'
    assert_allclose(torus.center, [0,0,0]) 

    torus.name = 'torus1'
    assert torus.name == 'torus1'
 
    torus.radius = 15
    assert torus.radius == 15

    torus.tube_radius = 4
    assert torus.tube_radius == 4
    
    torus.color = 'RED'
    assert torus.color == 'RED'

    torus.center = [4. 7. 9]
    assert_allclose(torus.center, [4, 7, 9])

    assert isinstance(torus, Shape)

def test_tube():
    point_list = [[2, 4, 5], [2, 6, 4], [1, 5, 8]]
    tube = Tube('tube', radius=10, points=point_list, \
                        color='RED', origin=[0,0,0])

    assert tube.name == 'tube' 			
    assert tube.radius == 10
    assert_allclose(tube.points, point_list)
    assert tube.color == 'RED'
    assert_allclose(tube.origin, [0,0,0]) 

    tube.name = 'tube1'
    assert tube.name == 'tube1'

    tube.radius = 15
    assert tube.radius == 15

    new_points_list = [[3, 4, 5], [1, 6 ,8], [2, 7, 3]
    tube.points = new_point_list
    assert_allclose(tube.points, new_point_list)
    
    tube.color = 'PINK'
    assert tube.color == 'PINK'
    
    tube.origin = [4. 7. 9]
    assert_allclose(tube.origin, [4, 7, 9])

    assert isinstance(tube, Shape)
def test_torus_knot():
    pass

 
