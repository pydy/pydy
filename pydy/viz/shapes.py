#!/usr/bin/env python

import numpy as np
try:
    import pythreejs as p3js
except ImportError:
    p3js = None

__all__ = [
           'Box',
           'Circle',
           'Cone',
           'Cube',
           'Cylinder',
           'Icosahedron',
           'Octahedron',
           'Plane',
           'Sphere',
           'Tetrahedron',
           'Torus',
           'TorusKnot',
           'Tube'
           ]


# This is a list of ColorKeywords from THREE.js
THREE_COLORKEYWORDS = ['aliceblue', 'antiquewhite', 'aqua',
                       'aquamarine', 'azure', 'beige', 'bisque',
                       'black', 'blanchedalmond', 'blue', 'blueviolet',
                       'brown', 'burlywood', 'cadetblue', 'chartreuse',
                       'chocolate', 'coral', 'cornflowerblue',
                       'cornsilk', 'crimson', 'cyan', 'darkblue',
                       'darkcyan', 'darkgoldenrod', 'darkgray',
                       'darkgreen', 'darkgrey', 'darkkhaki',
                       'darkmagenta', 'darkolivegreen', 'darkorange',
                       'darkorchid', 'darkred', 'darksalmon',
                       'darkseagreen', 'darkslateblue', 'darkslategray',
                       'darkslategrey', 'darkturquoise', 'darkviolet',
                       'deeppink', 'deepskyblue', 'dimgray', 'dimgrey',
                       'dodgerblue', 'firebrick', 'floralwhite',
                       'forestgreen', 'fuchsia', 'gainsboro',
                       'ghostwhite', 'gold', 'goldenrod', 'gray',
                       'green', 'greenyellow', 'grey', 'honeydew',
                       'hotpink', 'indianred', 'indigo', 'ivory',
                       'khaki', 'lavender', 'lavenderblush',
                       'lawngreen', 'lemonchiffon', 'lightblue',
                       'lightcoral', 'lightcyan',
                       'lightgoldenrodyellow', 'lightgray',
                       'lightgreen', 'lightgrey', 'lightpink',
                       'lightsalmon', 'lightseagreen', 'lightskyblue',
                       'lightslategray', 'lightslategrey',
                       'lightsteelblue', 'lightyellow', 'lime',
                       'limegreen', 'linen', 'magenta', 'maroon',
                       'mediumaquamarine', 'mediumblue',
                       'mediumorchid', 'mediumpurple', 'mediumseagreen',
                       'mediumslateblue', 'mediumspringgreen',
                       'mediumturquoise', 'mediumvioletred',
                       'midnightblue', 'mintcream', 'mistyrose',
                       'moccasin', 'navajowhite', 'navy', 'oldlace',
                       'olive', 'olivedrab', 'orange', 'orangered',
                       'orchid', 'palegoldenrod', 'palegreen',
                       'paleturquoise', 'palevioletred', 'papayawhip',
                       'peachpuff', 'peru', 'pink', 'plum',
                       'powderblue', 'purple', 'red', 'rosybrown',
                       'royalblue', 'saddlebrown', 'salmon',
                       'sandybrown', 'seagreen', 'seashell', 'sienna',
                       'silver', 'skyblue', 'slateblue', 'slategray',
                       'slategrey', 'snow', 'springgreen', 'steelblue',
                       'tan', 'teal', 'thistle', 'tomato', 'turquoise',
                       'violet', 'wheat', 'white', 'whitesmoke',
                       'yellow', 'yellowgreen']

MATERIALS = ["default", "checkerboard", "metal", "dirt", "foil", "water",
             "grass"]


class Shape(object):
    """Instantiates a shape. This is primarily used as a superclass for more
    specific shapes like Cube, Cylinder, Sphere etc.

    Shapes must be associated with a reference frame and a point using the
    VisualizationFrame class.

    Parameters
    ==========

    name : str, optional
        A name assigned to the shape.
    color : str, optional
        A color string from list of colors in THREE_COLORKEYWORDS

    Examples
    ========

    >>> from pydy.viz.shapes import Shape
    >>> s = Shape()
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> a = Shape(name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'

    """
    _p3js_material_attributes = {}

    def __init__(self, name='unnamed', color='grey', material="default"):

        self.name = name
        self.color = color
        self.material = material

        self.geometry_attrs = []

    def __str__(self):

        attributes = ([self.__class__.__name__,
                       self.name,
                       'color:' + self.color,
                       'material:' + self.material] +
                      sorted([attr + ':{}'.format(getattr(self, attr)) for
                              attr in self.geometry_attrs]))

        return ' '.join(['{}'] * len(attributes)).format(*attributes)

    def __repr__(self):
        return self.__class__.__name__

    @property
    def name(self):
        """Returns the name attribute of the shape."""
        return self._name

    @name.setter
    def name(self, new_name):
        """Sets the name attribute of the shape."""
        if not isinstance(new_name, str):
            raise TypeError("'name' should be a valid str object.")
        else:
            self._name = new_name

    @property
    def color(self):
        """Returns the color attribute of the shape."""
        return self._color

    @color.setter
    def color(self, new_color):
        """Sets the color attributes of the shape. This should be a valid
        three.js color keyword string."""
        if new_color not in THREE_COLORKEYWORDS:
            msg = "'color' should be a valid Three.js colors string:\n{}"
            raise ValueError(msg.format('\n'.join(THREE_COLORKEYWORDS)))
        else:
            self._color = new_color

    @property
    def material(self):
        """Returns the material attribute of the shape."""
        return self._material

    @material.setter
    def material(self, new_material):
        """Sets the material attribute of the shape, i.e. its shine,
        brightness, opacity etc.. The material should be a valid material
        from the listed MATERIALS. If a shape is attributed as "red" color,
        and "water" material, ideally it should have opacity and brightness
        properties like that of a red fluid.

        """
        if new_material.lower() not in MATERIALS:
            msg = "'material' is not valid. Choose from:\n{}"
            raise ValueError(msg.format('\n'.join(MATERIALS)))
        else:
            self._material = new_material

    def generate_dict(self, constant_map={}):
        """Returns a dictionary containing all the data associated with the
        Shape.

        Parameters
        ==========
        constant_map : dictionary
            If any of the shape's geometry are defined as SymPy expressions,
            then this dictionary should map all SymPy Symbol's found in the
            expressions to floats.

        """
        data_dict = {}
        data_dict['name'] = self.name
        data_dict['color'] = self.color
        data_dict['material'] = self.material
        data_dict['type'] = self.__repr__()
        for geom in self.geometry_attrs:
            atr = getattr(self, geom)
            try:
                data_dict[geom] = float(atr.subs(constant_map))
            except AttributeError:
                # not a SymPy expression
                data_dict[geom] = atr
            except TypeError:
                # can't convert expression to float
                raise TypeError('{} is an expression, you '.format(atr) +
                                'must provide a mapping to numerical values.')
        return data_dict

    def _p3js_mesh(self, constant_map={}):
        """Returns a PyThreeJS mesh object that corresponds to this shape."""

        if p3js is None:
            raise ImportError('pythreejs is not installed.')

        data = self.generate_dict(constant_map=constant_map)

        attrs = dict()

        for k, v in self._p3js_attribute_map.items():
            try:
                attrs[k] = data[v]
            except KeyError:
                attrs[k] = v

        try:
            geom_attr = '{}BufferGeometry'.format(self._p3js_geometry_type)
            Geometry = getattr(p3js, geom_attr)
        except AttributeError:
            geom_attr = '{}Geometry'.format(self._p3js_geometry_type)
            Geometry = getattr(p3js, geom_attr)

        geometry = Geometry(**attrs)

        # NOTE : For some reason traitlets doesn't think 'grey' is a valid HTML
        # color. This workaround should be removed, but deprecation will likely
        # be needed.
        if data['color'] == 'grey':
            color = 'gray'
        else:
            color = data['color']

        material = p3js.MeshStandardMaterial(color=color,
                                             **self._p3js_material_attributes)

        mesh = p3js.Mesh(name=data['name'], geometry=geometry,
                         material=material)

        self._mesh = mesh

        return mesh


class Box(Shape):
    """Instantiates a box of a given size.

    Parameters
    ==========
    width : float or SymPy expression
        Width of the box along the X axis.
    height : float or SymPy expression
        Height of the box along the Y axis.
    depth : float or SymPy expression
        Depth of the box along the Z axis.

    Examples
    ========

    >>> from pydy.viz.shapes import Box
    >>> s = Box(10.0, 5.0, 1.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.width
    5.0
    >>>s.height
    1.0
    >>>s.depth
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'

    """
    _p3js_geometry_type = 'Box'
    _p3js_attribute_map = {'width': 'width',
                           'height': 'height',
                           'depth': 'depth'}

    def __init__(self, width, height, depth, **kwargs):
        super(Box, self).__init__(**kwargs)
        self.geometry_attrs += ['width', 'height', 'depth']
        self.width = width
        self.height = height
        self.depth = depth


class Cube(Shape):
    """Instantiates a cube of a given size.

    Parameters
    ==========
    length : float or SymPy expression
        The length of the cube.

    Examples
    ========

    >>> from pydy.viz.shapes import Cube
    >>> s = Cube(10.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.length
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12.0
    >>> s.length
    12.0
    >>> a = Cube('my-shape2', 'red', length=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10.0

    """
    _p3js_geometry_type = 'Box'
    _p3js_attribute_map = {'width': 'length',
                           'height': 'length',
                           'depth': 'length'}

    def __init__(self, length, **kwargs):
        super(Cube, self).__init__(**kwargs)
        self.geometry_attrs.append('length')
        self.length = length


class Cylinder(Shape):
    """Instantiates a cylinder with given length and radius.

    Parameters
    ==========
    length : float or SymPy expression
        Length of the cylinder along its Y axis.
    radius : float or SymPy expression
        Radius of the cylinder (of the circular cross section normal to the Y
        axis).

    Examples
    ========

    >>> from pydy.viz.shapes import Cylinder
    >>> s = Cylinder(10.0, 5.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.length
    10.0
    >>> s.radius
    5.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12.0
    >>> s.length
    12.0
    >>> s.radius = 6.0
    >>> s.radius
    6.0
    >>> a = Cylinder(10.0, 5.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10.0
    >>> a.radius
    5.0

    """
    _p3js_geometry_type = 'Cylinder'
    _p3js_attribute_map = {'radiusTop': 'radius',
                           'radiusBottom': 'radius',
                           'radialSegments': 100,
                           'height': 'length'}

    def __init__(self, length, radius, **kwargs):
        super(Cylinder, self).__init__(**kwargs)
        self.geometry_attrs += ['length', 'radius']
        self.length = length
        self.radius = radius


class Cone(Shape):
    """Instantiates a cone with given length and base radius.

    Parameters
    ==========
    length : float or SymPy expression
        The length of the cone.
    radius : float or SymPy expression
        The base radius of the cone.

    Examples
    ========

    >>> from pydy.viz.shapes import Cone
    >>> s = Cone(10.0, 5.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.length
    10.0
    >>> s.radius
    5.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12.0
    >>> s.length
    12.0
    >>> s.radius = 6.0
    >>> s.radius
    6.0
    >>> a = Cone(10.0, 5.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10.0
    >>> a.radius
    5.0

    """
    _p3js_geometry_type = 'Cylinder'
    _p3js_attribute_map = {'radiusTop': 0.0,
                           'radiusBottom': 'radius',
                           'height': 'length'}

    def __init__(self, length, radius, **kwargs):
        super(Cone, self).__init__(**kwargs)
        self.geometry_attrs += ['length', 'radius']
        self.length = length
        self.radius = radius


class Sphere(Shape):
    """Instantiates a sphere with a given radius.

    Parameters
    ==========
    radius : float or SymPy expression
        The radius of the sphere.

    Examples
    ========

    >>> from pydy.viz.shapes import Sphere
    >>> s = Sphere(10.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.radius
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> a = Sphere(10.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0

    """
    _p3js_geometry_type = 'Sphere'
    _p3js_attribute_map = {'radius': 'radius',
                           'widthSegments': 100,
                           'heightSegments': 100}

    def __init__(self, radius=10.0, **kwargs):
        super(Sphere, self).__init__(**kwargs)
        self.geometry_attrs += ['radius']
        self.radius = radius


class Circle(Sphere):
    """Instantiates a circle with a given radius.

    Parameters
    ==========
    radius : float or SymPy Expression
        The radius of the circle.

    Examples
    ========

    >>> from pydy.viz.shapes import Circle
    >>> s = Circle(10.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.radius
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> a = Circle(10.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0

    """
    _p3js_geometry_type = 'Circle'
    _p3js_attribute_map = {'radius': 'radius',
                           'segments': 100}


class Plane(Shape):
    """Instantiates a plane with a given length and width.

    Parameters
    ==========
    length : float or SymPy expression
        Length of the plane along the Y axis.
    width : float or SymPy expression
        Width of the plane along the X axis.

    Examples
    ========

    >>> from pydy.viz.shapes import Plane
    >>> s = Plane(10.0, 5.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.length
    10.0
    >>> s.width
    5.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12.0
    >>> s.length
    12.0
    >>> s.width = 6.0
    >>> s.width
    6.0
    >>> a = Plane(10.0, 5.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10.0
    >>> a.width
    5.0

    """
    _p3js_geometry_type = 'Plane'
    _p3js_attribute_map = {'width': 'width',
                           'height': 'length'}
    _p3js_material_attributes = {'side': 'DoubleSide'}

    def __init__(self, length=10.0, width=5.0, **kwargs):
        super(Plane, self).__init__(**kwargs)
        self.geometry_attrs += ['length', 'width']
        self.length = length
        self.width = width


class Tetrahedron(Sphere):
    """Instantiates a Tetrahedron inscribed in a given radius circle.

    Parameters
    ==========
    radius : float or SymPy expression
        The radius of the circum-scribing sphere of around the tetrahedron.

    Examples
    ========

    >>> from pydy.viz.shapes import Tetrahedron
    >>> s = Tetrahedron(10.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.radius
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> a = Tetrahedron(10.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0

    """
    _p3js_geometry_type = 'Tetrahedron'
    _p3js_attribute_map = {'radius': 'radius'}


class Octahedron(Sphere):
    """Instantiaties an Octahedron inscribed in a circle of the given
    radius.

    Parameters
    ==========
    radius : float or SymPy expression.
        The radius of the circum-scribing sphere around the octahedron.

    Examples
    ========

    >>> from pydy.viz.shapes import Octahedron
    >>> s = Octahedron(10.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.radius
    10.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> a = Octahedron(10.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0

    """
    _p3js_geometry_type = 'Octahedron'
    _p3js_attribute_map = {'radius': 'radius'}


class Icosahedron(Sphere):
    """Instantiates an icosahedron inscribed in a sphere of the given
    radius.

    Parameters
    ==========
    radius : float or a SymPy expression
        Radius of the circum-scribing sphere for Icosahedron

    Examples
    ========

    >>> from pydy.viz.shapes import Icosahedron
    >>> s = Icosahedron(10)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>>s.radius
    10.0
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12
    >>> a = Icosahedron(10.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0

    """
    _p3js_geometry_type = 'Icosahedron'
    _p3js_attribute_map = {'radius': 'radius'}


class Torus(Shape):
    """Instantiates a torus with a given radius and section radius.

    Parameters
    ==========
    radius : float or SymPy expression
        The radius of the torus.
    tube_radius : float or SymPy expression
        The radius of the torus tube.

    Examples
    ========

    >>> from pydy.viz.shapes import Torus
    >>> s = Torus(10.0, 5.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.radius
    10.0
    >>> s.tube_radius
    5.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> s.tube_radius = 6.0
    >>> s.tube_radius
    6.0
    >>> a = Torus(10.0, 5.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0
    >>> a.tube_radius
    5.0

    """
    _p3js_geometry_type = 'Torus'
    _p3js_attribute_map = {'radius': 'radius',
                           'tube': 'tube_radius'}

    def __init__(self, radius, tube_radius, **kwargs):
        super(Torus, self).__init__(**kwargs)
        self.geometry_attrs += ['radius', 'tube_radius']
        self.radius = radius
        self.tube_radius = tube_radius

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, new_radius):
        self._radius = new_radius

    @property
    def tube_radius(self):
        return self._tube_radius

    @tube_radius.setter
    def tube_radius(self, new_tube_radius):
        self._tube_radius = new_tube_radius


class TorusKnot(Torus):
    """Instantiates a torus knot with given radius and section radius.

    Parameters
    ==========
    radius : float or SymPy expression
        The radius of the torus knot.
    tube_radius : float or SymPy expression
        The radius of the torus knot tube.

    Examples
    ========

    >>> from pydy.viz.shapes import TorusKnot
    >>> s = TorusKnot(10.0, 5.0)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.radius
    10.0
    >>> s.tube_radius
    5.0
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12.0
    >>> s.radius
    12.0
    >>> s.tube_radius = 6.0
    >>> s.tube_radius
    6.0
    >>> a = TorusKnot(10.0, 5.0, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10.0
    >>> a.tube_radius
    5.0

    """
    _p3js_geometry_type = 'TorusKnot'
    _p3js_attribute_map = {'radius': 'radius',
                           'tube': 'tube_radius'}


class Tube(Shape):
    """Instantiates a tube that sweeps along a path.

    Parameters
    ==========
    radius : float or SymPy expression
        The radius of the tube.
    points : array_like, shape(n, 3)
        An array of n (x, y, z) coordinates representing points that the
        tube's center line should follow.

    Examples
    ========

    >>> from pydy.viz.shapes import Tube
    >>> points = [[1.0, 2.0, 1.0], [2.0, 1.0, 1.0], [2.0, 3.0, 4.0]]
    >>> s = Tube(10.0, points)
    >>> s.name
    'unnamed'
    >>> s.color
    'grey'
    >>> s.points
    [[1.0, 2.0, 1.0], [2.0, 1.0, 1.0], [2.0, 3.0, 4.0]]
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 14.0
    >>> s.radius
    14.0
    >>> s.points = [[2.0, 1.0, 4.0], [1.0, 2.0, 4.0],
    ...             [2.0, 3.0, 1.0], [1.0, 1.0, 3.0]]
    >>> s.points
    [[2.0, 1.0, 4.0], [1.0, 2.0, 4.0], [2.0, 3.0, 1.0], [1.0, 1.0, 3.0]]
    >>> a = Tube(12.0, points, name='my-shape2', color='red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    12.0
    >>> a.points
    [[1.0, 2.0, 1.0], [2.0, 1.0, 1.0], [2.0, 3.0, 4.0]]

    """

    def __init__(self, radius, points, **kwargs):
        super(Tube, self).__init__(**kwargs)
        self.geometry_attrs += ['radius', 'points']
        self.radius = radius
        self.points = points

    @property
    def points(self):
        return self._points.tolist()

    @points.setter
    def points(self, new_points):
        self._points = np.asarray(new_points)
