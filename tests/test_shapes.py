from sympy.physics.mechanics import *
from numpy.testing import assert_allclose
#from pydy-viz import Shape, Cube, 

def test_cube():
    cube = Cube('cube', length=10, color='BLUE')
    assert cube.name == 'cube'
    assert cube.__str__ == 'Cube BLUE 10'
    assert cube.__repr__ == 'Cube'
    assert cube.__repr__ == 'cube'
    assert cube.length == 10
    assert cube.color == 'BLUE'
    
    cube.name = 'cube1'
    assert cube.name == 'cube1'
    
    cube.length = 16
    assert cube.length == 16

    cube.color = 'RED'
    assert cube.color == 'RED'
    
    cube.generate_dict() == {"color": (1.0, 0.0, 0.0), "type": "Cube", \
                                   "name": "cube1", "length": 16}

    assert isinstance(cube, Shape)

def test_cylinder():
    cylinder = Cylinder('cylinder', height=10, radius=5, color='blue')
    assert cylinder.name == 'cylinder'
    assert cylinder.__str__ == 'Cylinder blue height=10,radius=5'
    assert cylinder.__repr__ == 'Cylinder'
    assert cylinder.height == 10
    assert cylinder.radius == 5
    assert cylinder.color == 'blue'
    
    cylinder.name = 'cylinder1'
    assert cylinder.name == 'cylinder1'
    
    cylinder.height = 14
    assert cylinder.height == 14

    cylinder.radius = 7
    assert cylinder.radius == 7

    cylinder.color = 'cyan'
    assert cylinder.color == 'cyan'

    cylinder.generate_dict() == {"color": (0.0, 1.0, 1.0), \
                                         "type": "Cylinder", \
                                   "name": "cylinder1", "height": 14, \
                                       "radius" : 7}
                                       
    assert isinstance(cylinder, Shape)
    
def test_cone():
    cone = Cone('cone',heigth=10, radius=5, color='darkblue')
    assert cone.name == 'cone'
    assert cone.__str__ == 'Cone height=10,radius=5'
    assert cone.__repr__ == 'Cone'
    assert cone.height == 10
    assert cone.radius == 5
    assert cone.color == 'darkblue'
    
    cone.name = 'cone1'
    assert cone.name == 'cone1'
    
    cone.height = 16
    assert cone.height == 16

    cone.radius = 3
    assert cone.radius == 3

    cone.color = 'darkcyan'
    assert cone.color == 'darkcyan'                                            
    
    cone.generate_dict() == {"color": (0.0, 0.5450980392156862, \
                                            0.5450980392156862), \
                                              "type": "Cone", \
                                       "name": "cone1", "height": 16, \
                                       "radius" : 3}
    assert isinstance(cone, Shape)
                            
def test_sphere():
    sphere = Sphere('sphere', radius=10, color='azure')
    assert sphere.name == 'sphere'
    assert sphere.__str__ == 'Sphere azure 10'
    assert sphere.__repr__ == 'Sphere'
    assert sphere.radius == 10
    assert sphere.color == 'azure'
    
    sphere.name = 'sphere1'
    assert sphere.name == 'sphere1'
    
    sphere.radius = 14
    assert sphere.radius == 14

    sphere.color = 'aqua'
    assert sphere.color == 'aqua'

    sphere.generate_dict() == {"color":  (0.0, 1.0, 1.0), \
                                            "type": "Sphere", \
                                      "name": "sphere1", "radius" : 14}
    assert isinstance(sphere, Shape)
 
def test_circle():
    circle = Circle('circle', radius=10, color='gold')

    assert circle.name == 'circle'
    assert circle.__str__ == 'Circle gold 10'
    assert circle.__repr__ == 'Circle'
    assert circle.radius == 10
    assert circle.color == 'gold'
    
    circle.name = 'circle1'
    assert circle.name == 'circle1'
    
    circle.radius = 12
    assert circle.radius == 12

    circle.color = 'black'
    assert circle.color == 'black'

    circle.generate_dict() == {"color":  (0.0, 0.0, 0.0), \
                                         "type": "Circle", \
                                   "name": "circle1", "radius" : 12}
                                   
    assert isinstance(circle, Shape)

def test_mesh():
    point_list = [[2., 3., 1.], [4., 6., 2.], \
                         [5., 3., 1.], [5., 3., 6.], \
                         [2., 8., 4.], [7., 4., 1.]]
    
    mesh_shape = Mesh('mesh', points=point_list, color='green')                     
    assert mesh_shape.name == 'mesh'
    assert mesh_shape.__str__ == 'Mesh green'
    assert mesh_shape.__repr__ == 'Mesh'
    assert_allclose(mesh_shape.points,point_list)
    assert mesh_shape.color == 'green'
    
    
    mesh_shape.name = 'mesh1'
    assert mesh_shape.name == 'mesh1'
    
    
    new_point_list = [[3., 4., 12.], [2., 4., 4.], [3., 2., 41.], [2., 5., 4.]]
    mesh_shape.points = new_point_list
    assert_allclose(mesh_shape.points, new_point_list)
    
    mesh_shape.color = 'pink'
    assert mesh_shape.color == 'pink'
    
    mesh_shape.generate_dict() == {"color":  (1.0, 0.7529411764705882, \
                                          0.796078431372549), \
                                               "type": "Mesh", \
                                        "name": "mesh_shape1", \
                                           "points" : new_points_list}
                                          
    assert isinstance(mesh_shape, Shape)                      

def test_plane():
    plane = Plane('plane', length=10, width=20, color='indigo')
    assert plane.name == 'plane'
    assert plane.__str__ == 'Plane indigo 10 * 20'
    assert plane.__repr__ == 'plane'
    assert plane.length == 10
    assert plane.width == 20
    assert plane.color == 'indigo'
    
    plane.name = 'plane1'
    assert plane.name == 'plane1'
    
    plane.length = 30
    assert plane.length == 30
    
    plane.width = 10
    assert plane.width == 10
    
    plane.color = 'lavender'
    assert plane.color == 'lavender'
    
    plane.generate_dict() == {"color":  (0.9019607843137255, \
                                    0.9019607843137255, \
                                    0.9803921568627451), \
                                         "type": "Plane", \
                          "name": "plane1", "width" : 10, "length" : 30}
                                   
    assert isinstance(plane, Shape)

def test_tetrahedron():
    #Tetrahedron,Octahedron and Icosahedron
    # geometry is defined by the radius of the
    #circumscribed sphere. It would be mentioned explicitly in the
    #docstrings
    tetrahedron = Tetrahedron('tetrahedron', radius=5, color='maroon')
    assert tetrahedron.name == 'tetrahedron'
    assert tetrahedron.__str__ == 'Tetrahedron maroon 5'
    assert tetrahedron.__repr__ == 'Tetrahedron'
    assert tetrahedron.radius == 5
    assert tetrahedron.color == 'maroon'
    
    tetrahedron.name = 'tetrahedron1'
    assert tetrahedron.name == 'tetrahedron1'
    
    tetrahedron.radius = 7
    assert tetrahedron.radius == 7

    tetrahedron.color = 'orange'
    assert tetrahedron.color == 'orange'

    tetrahedron.generate_dict() == {"color":  (1.0, \
                                        0.6470588235294118, 0.0), \
                                          "type": "Tetrahedron", \
                                   "name": "tetrahedron1", "radius" : 7}
    assert isinstance(tetrahedron, Shape)

def test_octahedron():
    octahedron = Octahedron('octahedron', radius=12, color='purple')
    assert octahedron.name == 'octahedron'
    assert octahedron.__str__ == 'Octahedron purple 12'
    assert octahedron.__repr__ == 'Octahedron'
    assert octahedron.radius == 12
    assert octahedron.color == 'purple'
    
    octahedron.name = 'octahedron1'
    assert octahedron.name == 'octahedron1'
    
    octahedron.radius = 2
    assert octahedron.radius == 2

    octahedron.color = 'red'
    assert octahedron.color == 'red'

    octahedron.generate_dict() == {"color": (1.0, 0.0, 0.0), \
                                          "type": "Octahedron", \
                                   "name": "octahedron1", "radius" : 2}
                                   
    assert isinstance(octahedron, Shape)

def test_icosahedron():
    icosahedron = Icosahedron('icosahedron', radius=11, color='#FDF5E6')
    assert icosahedron.name == 'icosahedron #FDF5E6 11'
    assert icosahedron.__str__ == 'Icosahedron '
    assert icosahedron.__repr__ == 'Icosahedron'
    assert icosahedron.radius == 11
    assert icosahedron.color == '#FDF5E6'
    
    icosahedron.name = 'icosahedron1'
    assert icosahedron.name == 'icosahedron1'
    
    icosahedron.radius = 3
    assert icosahedron.radius == 3

    icosahedron.color = '#FFC0CB'
    assert icosahedron.color == '#FFC0CB'

    icosahedron.generate_dict() == {"color": (1.0, \
                                    0.7529411764705882, \
                                          0.796078431372549), \
                                          "type": "Icosahedron", \
                                   "name": "icosahedron1", "radius" : 3}
                                   
    assert isinstance(icosahedron, Shape)

def test_torus():
    torus = Torus('torus', radius=10, tube_radius=2, color='#FFFF00')

    assert torus.name == 'torus' 			
    assert torus.__str__ == 'Torus #FFFF00 radius=10,tube radius=2' 			
    assert torus.__repr__ == 'Torus' 			
    assert torus.radius == 10
    assert torus.tube_radius == 2
    assert torus.color == '#FFFF00'
    
    torus.name = 'torus1'
    assert torus.name == 'torus1'
    
    torus.radius = 15
    assert torus.radius == 15

    torus.tube_radius = 4
    assert torus.tube_radius == 4
    
    torus.color = '#FFFFFF'
    assert torus.color == '#FFFFFF'

    torus.generate_dict() == {"color": (1.0, 1.0, 1.0), \
                                          "type": "Torus", \
                                   "name": "torus1", "radius" : 15, \
                                        "tube_radius" : 4}
                                        
    assert isinstance(torus, Shape)

def test_tube():
    point_list = [[2., 4., 5.], [2., 6., 4.], [1., 5., 8.]]
    tube = Tube('tube', radius=10, points=point_list, color='#4682B4')

    assert tube.name == 'tube' 			
    assert tube.__str__ == 'Tube #4682B4 10' 			
    assert tube.__repr__ == 'Tube' 			
    assert tube.radius == 10
    assert_allclose(tube.points, point_list)
    assert tube.color == '#4682B4'
    

    tube.name = 'tube1'
    assert tube.name == 'tube1'
    	

    tube.radius = 15
    assert tube.radius == 15

    new_points_list = [[3., 4., 5.], [1, 6., 8.], [2., 7., 3.]]
    tube.points = new_point_list
    assert_allclose(tube.points, new_point_list)
    
    tube.color = 'pink'
    assert tube.color == 'pink'
    
    tube.generate_dict() == {"color": (1.0, \
                                      0.7529411764705882, \
                                        0.796078431372549), \
                                          "type": "Tube", \
                                   "name": "tube1", "radius" : 3, \
                                    "points" : new_points_list}
    
    assert isinstance(tube, Shape)
    
def test_torus_knot():
    
    torus_knot = TorusKnot('torus_knot', radius=10, tube_radius=2, color='#C0C0C0')

    assert torus_knot.name == 'torus_knot' 			
    assert torus_knot.__str__ == 'TorusKnot #C0C0C0 radius=10,tube radius=2' 			
    assert torus_knot.__repr__ == 'TorusKnot' 			
    assert torus_knot.radius == 10
    assert torus_knot.tube_radius == 2
    assert torus_knot.color == '#C0C0C0'
    
    torus_knot.name = 'torus_knot1'
    assert torus_knot.name == 'torus_knot1'
    	
    torus_knot.radius = 12
    assert torus_knot.radius == 12

    torus_knot.tube_radius = 1
    assert torus_knot.tube_radius == 1
    
    torus_knot.color = ''
    assert torus_knot.color == '#2E8B57'

    torus_knot.generate_dict() == {"color": (0.1803921568627451, \
                                              0.5450980392156862, \
                                              0.3411764705882353), \
                                          "type": "TorusKnot", \
                                   "name": "torus_knot1", "radius" : 12, \
                                    "tube_radius" : 1}

    assert isinstance(torus_knot, Shape)
