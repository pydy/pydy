#!/usr/bin/env python

# external
import numpy as np
from sympy import symbols
from numpy.testing import assert_allclose

# local
from ..shapes import *


def test_shape():
    shape = Shape(name='shape', color='BLUE')

    assert shape.name == 'shape'
    assert shape.__str__() == 'Shape shape color:BLUE'
    assert shape.__repr__() == 'Shape'
    assert shape.color == 'BLUE'

    shape.name = 'shape1'
    assert shape.name == 'shape1'

    shape.color = 'RED'
    assert shape.color == 'RED'

    assert shape.generate_dict() == {"color": "RED",
                                     "type": "Shape",
                                     "name": "shape1"}

    assert isinstance(shape, Shape)

    #testing unnamed
    shape_ = Shape(color='BLUE')

    assert shape_.name == 'unnamed'
    assert shape_.__str__() == 'Shape unnamed color:BLUE'
    assert shape_.__repr__() == 'Shape'


def test_shape_geometry_with_expressions():

    shape = Shape()
    shape.length = symbols('l')
    shape.geometry_attrs.append('length')
    expected = {"color": "grey",
                "type": "Shape",
                "name": "unnamed",
                "length": 16.0}
    actual = shape.generate_dict(constant_map={symbols('l'): 16.0})
    assert actual == expected

    shape = Shape()
    shape.length = symbols('l1') + symbols('l2') ** 2
    shape.geometry_attrs.append('length')
    expected = {"color": "grey",
                "type": "Shape",
                "name": "unnamed",
                "length": 20.0}
    actual = shape.generate_dict(constant_map={symbols('l1'): 4.0,
                                               symbols('l2'): 4.0})
    assert actual == expected


def test_cube():
    cube = Cube(10.0, name='cube', color='BLUE')

    assert cube.name == 'cube'
    assert cube.__str__() == 'Cube cube color:BLUE length:10.0'
    assert cube.__repr__() == 'Cube'
    assert cube.length == 10.0
    assert cube.color == 'BLUE'

    cube.name = 'cube1'
    assert cube.name == 'cube1'

    cube.length = 16.0
    assert cube.length == 16.0

    cube.color = 'RED'
    assert cube.color == 'RED'

    assert cube.generate_dict() == {"color": "RED",
                                    "type": "Cube",
                                    "name": "cube1",
                                    "length": 16.0}

    assert isinstance(cube, Shape)

    #testing unnamed
    cube = Cube(10.0, color='BLUE')

    assert cube.name == 'unnamed'
    assert cube.__str__() == 'Cube unnamed color:BLUE length:10.0'
    assert cube.__repr__() == 'Cube'

    cube = Cube(symbols('V') ** (1.0 / 3.0))
    actual = cube.generate_dict(constant_map={symbols('V'): 27.0})
    assert actual == {"color": "grey",
                      "type": "Cube",
                      "name": "unnamed",
                      "length": 3.0}


def test_cylinder():
    cylinder = Cylinder(10.0, 5.0, name='cylinder', color='blue')
    assert cylinder.name == 'cylinder'
    assert cylinder.__str__() == \
        'Cylinder cylinder color:blue length:10.0 radius:5.0'
    assert cylinder.__repr__() == 'Cylinder'
    assert cylinder.length == 10.0
    assert cylinder.radius == 5.0
    assert cylinder.color == 'blue'

    cylinder.name = 'cylinder1'
    assert cylinder.name == 'cylinder1'

    cylinder.length = 14.0
    assert cylinder.length == 14.0

    cylinder.radius = 7.0
    assert cylinder.radius == 7.0

    cylinder.color = 'cyan'
    assert cylinder.color == 'cyan'

    assert cylinder.generate_dict() == {"color": "cyan",
                                        "type": "Cylinder",
                                        "name": "cylinder1",
                                        "length": 14.0,
                                        "radius": 7.0}

    assert isinstance(cylinder, Shape)

    cylinder_ = Cylinder(10.0, 5.0, color='blue')
    assert cylinder_.name == 'unnamed'
    assert cylinder_.__str__() == \
        'Cylinder unnamed color:blue length:10.0 radius:5.0'

    assert cylinder_.__repr__() == 'Cylinder'


def test_cone():
    cone = Cone(10.0, 5.0, name='cone', color='darkblue')
    assert cone.name == 'cone'
    assert cone.__str__() == \
        'Cone cone color:darkblue length:10.0 radius:5.0'
    assert cone.__repr__() == 'Cone'
    assert cone.length == 10.0
    assert cone.radius == 5.0
    assert cone.color == 'darkblue'

    cone.name = 'cone1'
    assert cone.name == 'cone1'

    cone.length = 16.0
    assert cone.length == 16.0

    cone.radius = 3.0
    assert cone.radius == 3.0

    cone.color = 'darkcyan'
    assert cone.color == 'darkcyan'

    assert cone.generate_dict() == {"color": "darkcyan",
                                    "type": "Cone",
                                    "name": "cone1",
                                    "length": 16.0,
                                    "radius": 3.0}
    assert isinstance(cone, Shape)

    cone_ = Cone(10.0, 5.0, color='darkblue')
    assert cone_.name == 'unnamed'
    assert cone_.__str__() == \
        'Cone unnamed color:darkblue length:10.0 radius:5.0'
    assert cone_.__repr__() == 'Cone'


def test_sphere():
    sphere = Sphere(10.0, name='sphere', color='azure')
    assert sphere.name == 'sphere'
    assert sphere.__str__() == 'Sphere sphere color:azure radius:10.0'
    assert sphere.__repr__() == 'Sphere'
    assert sphere.radius == 10.0
    assert sphere.color == 'azure'

    sphere.name = 'sphere1'
    assert sphere.name == 'sphere1'

    sphere.radius = 14.0
    assert sphere.radius == 14.0

    sphere.color = 'aqua'
    assert sphere.color == 'aqua'

    assert sphere.generate_dict() == {"color":  "aqua",
                                      "type": "Sphere",
                                      "name": "sphere1",
                                      "radius": 14.0}
    assert isinstance(sphere, Shape)

    sphere_ = Sphere(10.0, color='azure')
    assert sphere_.name == 'unnamed'
    assert sphere_.__str__() == 'Sphere unnamed color:azure radius:10.0'
    assert sphere_.__repr__() == 'Sphere'


def test_circle():
    circle = Circle(10.0, name='circle', color='gold')

    assert circle.name == 'circle'
    assert circle.__str__() == 'Circle circle color:gold radius:10.0'
    assert circle.__repr__() == 'Circle'
    assert circle.radius == 10.0
    assert circle.color == 'gold'

    circle.name = 'circle1'
    assert circle.name == 'circle1'

    circle.radius = 12.0
    assert circle.radius == 12.0

    circle.color = 'black'
    assert circle.color == 'black'

    assert circle.generate_dict() == {"color": "black",
                                      "type": "Circle",
                                      "name": "circle1",
                                      "radius": 12.0}

    assert isinstance(circle, Shape)

    circle = Circle(10.0, color='gold')
    assert circle.name == 'unnamed'
    assert circle.__str__() == 'Circle unnamed color:gold radius:10.0'
    assert circle.__repr__() == 'Circle'


def test_mesh():
    point_list = [[2., 3., 1.], [4., 6., 2.],
                  [5., 3., 1.], [5., 3., 6.],
                  [2., 8., 4.], [7., 4., 1.]]

    mesh_shape = Mesh(point_list, name='mesh', color='green')
    assert mesh_shape.name == 'mesh'
    assert mesh_shape.__str__() == \
        'Mesh mesh color:green points:[[ 2.  3.  1.]\n [ 4.  6.  2.]\n [ 5.  3.  1.]\n [ 5.  3.  6.]\n [ 2.  8.  4.]\n [ 7.  4.  1.]]'
    assert mesh_shape.__repr__() == 'Mesh'
    assert_allclose(mesh_shape.points, point_list)
    assert mesh_shape.color == 'green'

    mesh_shape.name = 'mesh1'
    assert mesh_shape.name == 'mesh1'

    new_point_list = [[3., 4., 12.],
                      [2., 4., 4.],
                      [3., 2., 41.],
                      [2., 5., 4.]]
    mesh_shape.points = new_point_list
    assert_allclose(mesh_shape.points, new_point_list)

    mesh_shape.color = 'pink'
    assert mesh_shape.color == 'pink'

    actual = mesh_shape.generate_dict()

    expected = {"color": "pink",
                "type": "Mesh",
                "name": "mesh1",
                "points": np.asarray(new_point_list)}

    for key in ['color', 'type', 'name']:
        actual[key] == expected[key]
    assert_allclose(actual['points'], expected['points'])

    assert isinstance(mesh_shape, Shape)

    mesh_shape_ = Mesh(points=point_list, color='green')
    assert mesh_shape_.name == 'unnamed'
    assert mesh_shape_.__str__() == \
        'Mesh unnamed color:green points:[[ 2.  3.  1.]\n [ 4.  6.  2.]\n [ 5.  3.  1.]\n [ 5.  3.  6.]\n [ 2.  8.  4.]\n [ 7.  4.  1.]]'
    assert mesh_shape_.__repr__() == 'Mesh'


def test_plane():
    plane = Plane(10.0, 20.0, name='plane', color='indigo')
    assert plane.name == 'plane'
    assert plane.__str__() == \
        'Plane plane color:indigo length:10.0 width:20.0'
    assert plane.__repr__() == 'Plane'
    assert plane.length == 10.0
    assert plane.width == 20.0
    assert plane.color == 'indigo'

    plane.name = 'plane1'
    assert plane.name == 'plane1'

    plane.length = 30.0
    assert plane.length == 30.0

    plane.width = 10.0
    assert plane.width == 10.0

    plane.color = 'lavender'
    assert plane.color == 'lavender'

    assert plane.generate_dict() == {"color":  "lavender",
                                     "type": "Plane",
                                     "name": "plane1",
                                     "width": 10.0,
                                     "length": 30.0}

    assert isinstance(plane, Shape)

    plane_ = Plane(10.0, 20.0, color='indigo')
    assert plane_.name == 'unnamed'
    assert plane_.__str__() == \
        'Plane unnamed color:indigo length:10.0 width:20.0'
    assert plane_.__repr__() == 'Plane'


def test_tetrahedron():
    #Tetrahedron,Octahedron and Icosahedron
    # geometry is defined by the radius of the
    #circumscribed sphere. It would be mentioned explicitly in the
    #docstrings
    tetrahedron = Tetrahedron(5.0, name='tetrahedron', color='maroon')
    assert tetrahedron.name == 'tetrahedron'
    assert tetrahedron.__str__() == \
        'Tetrahedron tetrahedron color:maroon radius:5.0'
    assert tetrahedron.__repr__() == 'Tetrahedron'
    assert tetrahedron.radius == 5.0
    assert tetrahedron.color == 'maroon'

    tetrahedron.name = 'tetrahedron1'
    assert tetrahedron.name == 'tetrahedron1'

    tetrahedron.radius = 7.0
    assert tetrahedron.radius == 7.0

    tetrahedron.color = 'orange'
    assert tetrahedron.color == 'orange'

    assert tetrahedron.generate_dict() == {"color":  "orange",
                                           "type": "Tetrahedron",
                                           "name": "tetrahedron1",
                                           "radius": 7.0}
    assert isinstance(tetrahedron, Shape)

    tetrahedron_ = Tetrahedron(5.0, color='maroon')
    assert tetrahedron_.name == 'unnamed'
    assert tetrahedron_.__str__() == \
        'Tetrahedron unnamed color:maroon radius:5.0'
    assert tetrahedron_.__repr__() == 'Tetrahedron'


def test_octahedron():
    octahedron = Octahedron(12.0, name='octahedron', color='purple')
    assert octahedron.name == 'octahedron'
    assert octahedron.__str__() == \
        'Octahedron octahedron color:purple radius:12.0'
    assert octahedron.__repr__() == 'Octahedron'
    assert octahedron.radius == 12.0
    assert octahedron.color == 'purple'

    octahedron.name = 'octahedron1'
    assert octahedron.name == 'octahedron1'

    octahedron.radius = 2.0
    assert octahedron.radius == 2.0

    octahedron.color = 'red'
    assert octahedron.color == 'red'

    assert octahedron.generate_dict() == {"color": "red",
                                          "type": "Octahedron",
                                          "name": "octahedron1",
                                          "radius": 2.0}

    assert isinstance(octahedron, Shape)

    octahedron_ = Octahedron(12.0, color='purple')
    assert octahedron_.name == 'unnamed'
    assert octahedron_.__str__() == \
        'Octahedron unnamed color:purple radius:12.0'
    assert octahedron_.__repr__() == 'Octahedron'


def test_icosahedron():
    icosahedron = Icosahedron(11.0, name='icosahedron', color='#FDF5E6')
    assert icosahedron.name == 'icosahedron'
    assert icosahedron.__str__() == \
        'Icosahedron icosahedron color:#FDF5E6 radius:11.0'
    assert icosahedron.__repr__() == 'Icosahedron'
    assert icosahedron.radius == 11.0
    assert icosahedron.color == '#FDF5E6'

    icosahedron.name = 'icosahedron1'
    assert icosahedron.name == 'icosahedron1'

    icosahedron.radius = 3.0
    assert icosahedron.radius == 3.0

    icosahedron.color = '#FFC0CB'
    assert icosahedron.color == '#FFC0CB'

    assert icosahedron.generate_dict() == {"color": "#FFC0CB",
                                           "type": "Icosahedron",
                                           "name": "icosahedron1",
                                           "radius": 3.0}

    assert isinstance(icosahedron, Shape)

    icosahedron_ = Icosahedron(11.0, color='#FDF5E6')
    assert icosahedron_.name == 'unnamed'
    assert icosahedron_.__str__() == \
        'Icosahedron unnamed color:#FDF5E6 radius:11.0'
    assert icosahedron_.__repr__() == 'Icosahedron'


def test_torus():
    torus = Torus(10.0, 2.0, name='torus', color='#FFFF00')

    assert torus.name == 'torus'
    assert torus.__str__() == \
        'Torus torus color:#FFFF00 radius:10.0 tube_radius:2.0'
    assert torus.__repr__() == 'Torus'
    assert torus.radius == 10.0
    assert torus.tube_radius == 2.0
    assert torus.color == '#FFFF00'

    torus.name = 'torus1'
    assert torus.name == 'torus1'

    torus.radius = 15.0
    assert torus.radius == 15.0

    torus.tube_radius = 4.0
    assert torus.tube_radius == 4.0

    torus.color = '#FFFFFF'
    assert torus.color == '#FFFFFF'

    assert torus.generate_dict() == {"color": "#FFFFFF",
                                     "type": "Torus",
                                     "name": "torus1",
                                     "radius": 15.0,
                                     "tube_radius": 4.0}

    assert isinstance(torus, Shape)

    torus_ = Torus(10.0, 2.0, color='#FFFF00')
    assert torus_.name == 'unnamed'
    assert torus_.__str__() == \
        'Torus unnamed color:#FFFF00 radius:10.0 tube_radius:2.0'
    assert torus_.__repr__() == 'Torus'


def test_tube():
    point_list = [[2., 4., 5.], [2., 6., 4.], [1., 5., 8.]]
    tube = Tube(10.0, point_list, name='tube', color='#4682B4')

    assert tube.name == 'tube'
    assert tube.__str__() == 'Tube tube color:#4682B4 points:[[ 2.  4.  5.]\n [ 2.  6.  4.]\n [ 1.  5.  8.]] radius:10.0'
    assert tube.__repr__() == 'Tube'
    assert tube.radius == 10.0
    assert_allclose(tube.points, point_list)
    assert tube.color == '#4682B4'

    tube.name = 'tube1'
    assert tube.name == 'tube1'

    tube.radius = 15.0
    assert tube.radius == 15.0

    new_point_list = [[3., 4., 5.], [1, 6., 8.], [2., 7., 3.]]
    tube.points = new_point_list
    assert_allclose(tube.points, new_point_list)

    tube.color = 'pink'
    assert tube.color == 'pink'

    actual = tube.generate_dict()

    expected = {"color": "pink",
                "type": "Tube",
                "name": "tube1",
                "radius": 15.0,
                "points": np.asarray(new_point_list)}

    for key in ['color', 'type', 'name', 'radius']:
        actual[key] == expected[key]
    assert_allclose(actual['points'], expected['points'])

    assert isinstance(tube, Shape)

    tube_ = Tube(10.0, point_list, color='#4682B4')
    assert tube_.name == 'unnamed'
    assert tube_.__str__() == 'Tube unnamed color:#4682B4 points:[[ 2.  4.  5.]\n [ 2.  6.  4.]\n [ 1.  5.  8.]] radius:10.0'
    assert tube_.__repr__() == 'Tube'


def test_torus_knot():

    torus_knot = TorusKnot(10.0, 2.0, name='torus_knot', color='#C0C0C0')

    assert torus_knot.name == 'torus_knot'
    assert torus_knot.__str__() == \
        'TorusKnot torus_knot color:#C0C0C0 radius:10.0 tube_radius:2.0'
    assert torus_knot.__repr__() == 'TorusKnot'
    assert torus_knot.radius == 10.0
    assert torus_knot.tube_radius == 2.0
    assert torus_knot.color == '#C0C0C0'

    torus_knot.name = 'torus_knot1'
    assert torus_knot.name == 'torus_knot1'

    torus_knot.radius = 12.0
    assert torus_knot.radius == 12.0

    torus_knot.tube_radius = 1.0
    assert torus_knot.tube_radius == 1.0

    torus_knot.color = '#2E8B57'
    assert torus_knot.color == '#2E8B57'

    assert torus_knot.generate_dict() == {"color": "#2E8B57",
                                          "type": "TorusKnot",
                                          "name": "torus_knot1",
                                          "radius": 12.0,
                                          "tube_radius": 1}

    assert isinstance(torus_knot, Shape)

    torus_knot_ = TorusKnot(10.0, 2.0, color='#C0C0C0')
    assert torus_knot_.name == 'unnamed'
    assert torus_knot_.__str__() == \
        'TorusKnot unnamed color:#C0C0C0 radius:10.0 tube_radius:2.0'
    assert torus_knot_.__repr__() == 'TorusKnot'
