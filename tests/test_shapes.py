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
                                    color='blue', center=[0,0,0])
    assert cylinder.name == 'cylinder'
    assert cylinder.height == 10
    assert cylinder.radius == 5
    assert cylinder.color == 'blue'
    assert_allclose(cylinder.center, [0, 0, 0])
    
    
    cylinder.name = 'cylinder1'
    assert cylinder.name == 'cylinder1'

    cylinder.height = 14
    assert cylinder.height == 14

    cylinder.radius = 7
    assert cylinder.radius == 7

    cylinder.color = 'cyan'
    assert cylinder.color == 'cyan'

    cylinder.center = [3, 6, -2]
    assert_allclose(cylinder.center, [3, 6, -2])
    
    assert isinstance(cylinder, Shape)
    
def test_cone():
    cone = Cone('cone',heigth=10, radius=5, \
                        color='darkblue', center=[0,0,0])
    assert cone.name == 'cone'
    assert cone.height == 10
    assert cone.radius == 5
    assert cone.color == 'darkblue'
    assert_allclose(cone.center, [0, 0, 0])
    
    cone.name = 'cone1'
    assert cone.name == 'cone1'

    cone.height = 16
    assert cone.height == 16

    cone.radius = 3
    assert cone.radius == 3

    cone.color = 'darkcyan'
    assert cone.color == 'darkcyan'                                            
    
    cone.center = [3, 7, 3]
    assert_allclose(cone.center, [3, 7, 3])
    
    assert isinstance(cone, Shape)
                            
def test_sphere():
    sphere = Sphere('sphere', radius=10, color='azure', \
                                   center=[0,0,0])
    assert sphere.name == 'sphere1'
    assert sphere.radius == 10
    assert sphere.color == 'azure'
    assert_allclose(sphere.center, [0, 0, 0])	                                                       
    sphere.name = 'sphere1'
    assert sphere.name == 'sphere1'
   
    sphere.radius = 14
    assert sphere.radius == 14

    sphere.color = 'aqua'
    assert sphere.color == 'aqua'

    sphere.center = [5, 2, 6]
    assert_allclose(sphere.center, [5, 2, 6])
  
    assert isinstance(sphere, Shape)
 


def test_circle():
    circle = Circle('circle', radius=10, \
                           color='gold',center=[0,0,0])

    assert circle.name == 'circle'
    assert circle.radius == 10
    assert circle.color == 'gold'
    assert_allclose(circle.center, [0, 0, 0])

    circle.name = 'circle1'
    assert circle.name == 'circle1'

    circle.radius = 12
    assert circle.radius == 12

    circle.color = 'black'
    assert circle.color == 'black'

    circle.center = [0, -3, 5]
    assert_allclose(circle.center, [0, -3, 5])
    
    assert isinstance(circle, Shape)

def test_mesh_shape():
    point_list = [[2, 3, 1], [4, 6, 2], \
                         [5, 3, 1], [5, 3, 6], \
                         [2, 8, 4], [7, 4, 1]]
    
    mesh_shape = MeshShape('mesh_shape', points=point_list, \
                               color='green', origin=[0,0,0])                     
    assert mesh_shape.name == 'mesh_shape'
    assert_allclose(mesh_shape.points,point_list)
    assert mesh_shape.color == 'green'
    assert_allclose(mesh_shape.origin, point_list)
    
    mesh_shape.name = 'mesh_shape1'
    assert mesh_shape.name == 'mesh_shape1'
    
    new_point_list = [[3, 4, 12], [2, 4, 4], [3, 2, 41], [2, 5, 4]]
    mesh_shape.points = new_point_list
    assert_allclose(mesh_shape.points, new_point_list)
    
    mesh_shape.color = 'pink'
    assert mesh_shape.color == 'pink'
    
    mesh_shape.origin = [-5, 1, 2]
    assert_allclose(mesh_shape.origin, [-5, 1, 2])
    
    assert isinstance(mesh_shape, Shape)                      

def test_plane():
    plane = Plane('plane', length=10, width=20, \
                       color='indigo', center=[0,0,0])
    assert plane.name == 'plane'
    assert plane.length == 10
    assert plane.width == 20
    assert plane.color == 'indigo'
    assert_allclose(plane.center, [0, 0, 0])
    
    plane.name = 'plane1'
    assert plane.name == 'plane1'
    
    plane.length = 30
    assert plane.length == 30
    
    plane.width = 10
    assert plane.width == 10
    
    plane.color = 'lavender'
    assert plane.color == 'lavender'
    
    plane.center = [7, 4, 5]
    assert_allclose(plane.center, [7, 4, 5])
    
    assert isinstance(plane, Shape)

def test_tetrahedron():
    #Tetrahedron,Octahedron and Icosahedron
    # geometry is defined by the radius of the
    #circumscribed sphere. It would be mentioned explicitly in the
    #docstrings
    tetrahedron = Tetrahedron('tetrahedron', radius=5, \
                                    color='maroon', center=[0,0,0])
    assert tetrahedron.name == 'tetrahedron'
    assert tetrahedron.radius == 5
    assert tetrahedron.color == 'maroon'
    assert_allclose(tetrahedron.center, [0, 0, 0])
    
    
    tetrahedron.name = 'tetrahedron1'
    assert tetrahedron.name == 'tetrahedron1'

    tetrahedron.radius = 7
    assert tetrahedron.radius == 7

    tetrahedron.color = 'orange'
    assert tetrahedron.color == 'orange'

    tetrahedron.center = [0, 1, -4]
    assert_allclose(tetrahedron.center, [0, 1, -4])
    
    assert isinstance(tetrahedron, Shape)

def test_octahedron():
    octahedron = Octahedron('octahedron', radius=12, \
                                    color='purple', center=[0,0,0])
    assert octahedron.name == 'octahedron'
    assert octahedron.radius == 12
    assert octahedron.color == 'purple'
    assert_allclose(octahedron.center, [0, 0, 0])
    
    
    octahedron.name = 'octahedron1'
    assert octahedron.name == 'octahedron1'

    octahedron.radius = 2
    assert octahedron.radius == 2

    octahedron.color = 'red'
    assert octahedron.color == 'red'

    octahedron.center = [21, 11, 21]
    assert_allclose(octahedron.center, [21, 11, 21])
    
    assert isinstance(octahedron, Shape)

def test_icosahedron():
    icosahedron = Icosahedron('icosahedron', radius=11, \
                                    color='#FDF5E6', center=[0,0,0])
    assert icosahedron.name == 'icosahedron'
    assert icosahedron.radius == 11
    assert icosahedron.color == '#FDF5E6'
    assert_allclose(icosahedron.center, [0, 0, 0])
    
    
    icosahedron.name = 'icosahedron1'
    assert icosahedron.name == 'icosahedron1'

    icosahedron.radius = 3
    assert icosahedron.radius == 3

    icosahedron.color = '#FFC0CB'
    assert icosahedron.color == '#FFC0CB'

    icosahedron.center = [10, 11, 12]
    assert_allclose(icosahedron.center, [10, 11, 12])
    
    assert isinstance(icosahedron, Shape)

def test_torus():
    torus = Torus('torus', radius=10, tube_radius=2, \
                        color='#FFFF00', center=[0,0,0])

    assert torus.name == 'torus' 			
    assert torus.radius == 10
    assert torus.tube_radius == 2
    assert torus.color == '#FFFF00'
    assert_allclose(torus.center, [0,0,0]) 

    torus.name = 'torus1'
    assert torus.name == 'torus1'
 
    torus.radius = 15
    assert torus.radius == 15

    torus.tube_radius = 4
    assert torus.tube_radius == 4
    
    torus.color = '#FFFFFF'
    assert torus.color == '#FFFFFF'

    torus.center = [4. 7. 9]
    assert_allclose(torus.center, [4, 7, 9])

    assert isinstance(torus, Shape)

def test_tube():
    point_list = [[2, 4, 5], [2, 6, 4], [1, 5, 8]]
    tube = Tube('tube', radius=10, points=point_list, \
                        color='#4682B4', origin=[0,0,0])

    assert tube.name == 'tube' 			
    assert tube.radius == 10
    assert_allclose(tube.points, point_list)
    assert tube.color == '#4682B4'
    assert_allclose(tube.origin, [0,0,0]) 

    tube.name = 'tube1'
    assert tube.name == 'tube1'

    tube.radius = 15
    assert tube.radius == 15

    new_points_list = [[3, 4, 5], [1, 6 ,8], [2, 7, 3]
    tube.points = new_point_list
    assert_allclose(tube.points, new_point_list)
    
    tube.color = 'pink'
    assert tube.color == 'pink'
    
    tube.origin = [4. 7. 9]
    assert_allclose(tube.origin, [4, 7, 9])

    assert isinstance(tube, Shape)
    
def test_torus_knot():
    
    torus_knot = TorusKnot('torus_knot', radius=10, tube_radius=2, \
                        color='#C0C0C0', center=[0,0,0])

    assert torus_knot.name == 'torus_knot' 			
    assert torus_knot.radius == 10
    assert torus_knot.tube_radius == 2
    assert torus_knot.color == '#C0C0C0'
    assert_allclose(torus_knot.center, [0,0,0]) 

    torus_knot.name = 'torus_knot1'
    assert torus_knot.name == 'torus_knot1'
 
    torus_knot.radius = 12
    assert torus_knot.radius == 12

    torus_knot.tube_radius = 1
    assert torus_knot.tube_radius == 1
    
    torus_knot.color = '#2E8B57'
    assert torus_knot.color == '#2E8B57'

    torus_knot.center = [1. 2. 1]
    assert_allclose(torus_knot.center, [1, 2, 1])

    assert isinstance(torus_knot, Shape)
