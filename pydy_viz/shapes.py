#importing matplotlib's colorconverter here,
#module can be copied later in our source.
__all__ = ['Shape', \
           'Cube',  \
           'Cylinder', \
           'Cone', \
           'Sphere', \
           'Circle', \
           'Plane', \
           ]


from matplotlib.colors import ColorConverter
from sympy.physics.mechanics import Point, ReferenceFrame


convert = ColorConverter()

class Shape(object):
    """
    A Shape. It is a superclass for more general shapes like
    Mesh, Cylinder, Sphere etc. 
    
    Default Shape can be used for Particle visualizations, 
    as in sympy.physics.mechanics.particle
    
    The Shape classes are used for creating visualization objects, used
    for studying multibody dynamics and in animations.
    
    Values need to be supplied on initialization, but can be changed 
    later.

    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for this shape

    Examples
    ========

    >>> from pydy_viz.shapes import Shape
    >>> 
    >>> s = Shape()
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> a = Shape('my-shape2', 'red')
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey'):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        
    @property
    def name(self):
        """
        Name property of the shape,
        defines a name to a shape.
        """
        return self._name        
    
    @name.setter       
    def name(self, new_name):
        if not isinstance(new_name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = new_name    
    
    @property    
    def color(self):
        """    
        color property of the shape,
        used for visualizations
        """
        return self._color
  
    @color.setter
    def color(self, new_color):
        if not isinstance(new_color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = new_color
            self._color_rgb = convert.to_rgb(new_color)                             

    def color_in_rgb(self):
        """Returns the rgb value of the
           defined shape color.
        """
        return self._color_rgb

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        return self._data_dict
        
    def __str__(self):
        return 'Shape ' + self._name + ' color:' + self._color
    
    def __repr__(self):
        return 'Shape'    
        
class Cube(Shape):
    """
    A Cube. This class generates a Cube, with given length of side, 
    and color.Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Cube
    length: int or float.
        Length of side of Cube    

    Examples
    ========

    >>> from pydy_viz.shapes import Cube
    >>> 
    >>> s = Cube(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.length
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length
    12
    >>> a = Cube('my-shape2', 'red', length=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', length=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = length    
        
    def __str__(self):
        return 'Cube ' + self._name + ' color:' + self._color + \
                                        ' length:' + str(self._length)
    
    def __repr__(self):
        return 'Cube'    

    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = new_length

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cube, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['length'] = self._length
        
        return self._data_dict
        
        
class Cylinder(Shape):
    """
    A Cylinder. This class generates a Cylinder with given length,
    radius, and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to Cylinder
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Cylinder
    length: int or float, length of the Cylinder
    radius: int or float, radius of the Cylinder
        
    Examples
    ========

    >>> from pydy_viz.shapes import Cylinder
    >>> 
    >>> s = Cylinder(length=10, radius=5)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>> s.length
    10
    >>> s.radius
    5
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length
    12
    >>> s.radius = 6
    >>> s.radius
    6
    >>> a = Cylinder('my-shape2', 'red', length=10, radius=5)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10
    >>> a.radius
    5
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', \
                                  color='grey', length=10, radius=5):
                                      
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = length    
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius
            
    def __str__(self):
        return 'Cylinder ' + self._name + ' color:' + self._color + \
                                    ' length:' + str(self._length) + \
                                        ' radius:' + str(self._radius)
    def __repr__(self):
        return 'Cylinder'    

    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = new_length            
    
    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cylinder, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['length'] = self._length
        self._data_dict['radius'] = self._radius
        return self._data_dict

class Cone(Shape):
    """
    A Cone. This class generates a Cone with given length,
    base radius, and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to Cone
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Cone
    length: int or float, length of the Cone
    radius: int or float, base radius of the Cone
        
    Examples
    ========

    >>> from pydy_viz.shapes import Cone
    >>> 
    >>> s = Cone(length=10, radius=5)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>> s.length
    10
    >>> s.radius
    5
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length
    12
    >>> s.radius = 6
    >>> s.radius
    6
    >>> a = Cone('my-shape2', 'red', length=10, radius=5)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10
    >>> a.radius
    5
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', \
                                  color='grey', length=10, radius=5):
                                      
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = length    
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius
            
    def __str__(self):
        return 'Cone ' + self._name + ' color:' + self._color + \
                                    ' length:' + str(self._length) + \
                                        ' radius:' + str(self._radius)
    def __repr__(self):
        return 'Cone'    

    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = new_length            
    
    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cone, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['length'] = self._length
        self._data_dict['radius'] = self._radius
        return self._data_dict

class Sphere(Shape):
    """
    A Sphere. This class generates a Sphere, with given length of side, 
    and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Sphere
    radius: int or float.
        Radius of Sphere

    Examples
    ========

    >>> from pydy_viz.shapes import Sphere
    >>> 
    >>> s = Sphere(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.radius
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12
    >>> s.radius
    12
    >>> a = Sphere('my-shape2', 'red', radius=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', radius=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius    
        
    def __str__(self):
        return 'Sphere ' + self._name + ' color:' + self._color + \
                                        ' radius:' + str(self._radius)
    
    def __repr__(self):
        return 'Sphere'    

    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cube, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['radius'] = self._radius
        
        return self._data_dict

class Circle(Shape):
    """
    A Circle. This class generates a Circle, with given radius, 
    and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Circle
    radius: int or float.
        Radius of Circle

    Examples
    ========

    >>> from pydy_viz.shapes import Circle
    >>> 
    >>> s = Circle(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.radius
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12
    >>> s.radius
    12
    >>> a = Circle('my-shape2', 'red', radius=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', radius=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius    
        
    def __str__(self):
        return 'Circle ' + self._name + ' color:' + self._color + \
                                        ' radius:' + str(self._radius)
    
    def __repr__(self):
        return 'Circle'    

    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Circle, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['radius'] = self._radius
        
        return self._data_dict

class Plane(Shape):
    """
    A Plane. This class generates a Plane with given length,
    width, and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to Plane
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Plane
    length: int or float, length of the Plane
    width: int or float, radius of the Plane
        
    Examples
    ========

    >>> from pydy_viz.shapes import Plane
    >>> 
    >>> s = Plane(length=10, width=5)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>> s.length
    10
    >>> s.width
    5
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length
    12
    >>> s.width = 6
    >>> s.width
    6
    >>> a = Plane('my-shape2', 'red', length=10, width=5)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.length
    10
    >>> a.width
    5
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', \
                                  color='grey', length=10, width=5):
                                      
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = length    
        
        if not isinstance(width, (int, float)):
            raise TypeError('''Width should be a float or int''')
        else:
            self._width = width
            
    def __str__(self):
        return 'Plane ' + self._name + ' color:' + self._color + \
                                    ' length:' + str(self._length) + \
                                        ' width:' + str(self._width)
    def __repr__(self):
        return 'Plane'    

    @property
    def length(self):
        return self._length
    
    @length.setter
    def length(self, new_length):
        if not isinstance(new_length, (int, float)):
            raise TypeError('''Length should be a float or int''')
        else:
            self._length = new_length            
    
    @property
    def width(self):
        return self._width
    
    @width.setter
    def width(self, new_width):
        if not isinstance(new_width, (int, float)):
            raise TypeError('''Width should be a float or int''')
        else:
            self._width = new_width

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Plane, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['length'] = self._length
        self._data_dict['width'] = self._width
        return self._data_dict

class Tetrahedron(Shape):
    """
    A Tetrahedron. This class generates a Tetrahedron.
    The argument given is the radius of the circumscribing sphere of the 
    tetrahedron, 
    and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Tetrahedron
    radius: int or float.
        Radius of circum-scribing sphere of Tetrahedron

    Examples
    ========

    >>> from pydy_viz.shapes import Tetrahedron
    >>> 
    >>> s = Tetrahedron(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.radius
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12
    >>> s.radius
    12
    >>> a = Tetrahedron('my-shape2', 'red', radius=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', radius=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius    
        
    def __str__(self):
        return 'Tetrahedron ' + self._name + ' color:' + self._color + \
                                        ' radius:' + str(self._radius)
    
    def __repr__(self):
        return 'Tetrahedron'    

    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cube, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['radius'] = self._radius
        
        return self._data_dict

class Octahedron(Shape):
    """
    A Octahedron. This class generates a Octahedron.
    The argument given is the radius of the circumscribing sphere of the 
    octahedron, 
    and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Octahedron
    radius: int or float.
        Radius of the circum-scribing sphere for Octahedron

    Examples
    ========

    >>> from pydy_viz.shapes import Octahedron
    >>> 
    >>> s = Octahedron(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.radius
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12
    >>> s.radius
    12
    >>> a = Octahedron('my-shape2', 'red', radius=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', radius=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius    
        
    def __str__(self):
        return 'Octahedron ' + self._name + ' color:' + self._color + \
                                        ' radius:' + str(self._radius)
    
    def __repr__(self):
        return 'Octahedron'    

    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cube, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['radius'] = self._radius
        
        return self._data_dict

class Icosahedron(Shape):
    """
    A Icosahedron. This class generates a Icosahedron.
    The argument given is the radius of the circumscribing sphere of the 
    icosahedron, 
    and color. Default color is grey.
    
    Parameters
    ==========
    name : str
        Name assigned to shape
    color: str
        A color string from list of colors in pydy_viz.colors module
        This color is used in drawing visualizations for Icosahedron
    radius: int or float.
        Radius of the circum-scribing sphere for Icosahedron

    Examples
    ========

    >>> from pydy_viz.shapes import Icosahedron
    >>> 
    >>> s = Icosahedron(10)
    >>> s.name
    'UnNamed'
    >>> s.color
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.radius
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.radius = 12
    >>> s.radius
    12
    >>> a = Icosahedron('my-shape2', 'red', radius=10)
    >>> a.name
    'my-shape2'
    >>> a.color
    'red'
    >>> a.radius
    10
    >>> a.color_in_rgb()
    (1.0, 0.0, 0.0)
    
    """
    
    def __init__(self, name='UnNamed', color='grey', radius=10):
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object.')
        else:
            self._name = name
        
        if not isinstance(color, str):
            raise TypeError('''color should be a valid \
                               colors string. for info on colors, see \
                               pydy_viz.colors module''')       
        else:
            self._color = color
            self._color_rgb = convert.to_rgb(color)                       
        
        if not isinstance(radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = radius    
        
    def __str__(self):
        return 'Icosahedron ' + self._name + ' color:' + self._color + \
                                        ' radius:' + str(self._radius)
    
    def __repr__(self):
        return 'Icosahedron'    

    @property
    def radius(self):
        return self._radius
    
    @radius.setter
    def radius(self, new_radius):
        if not isinstance(new_radius, (int, float)):
            raise TypeError('''Radius should be a float or int''')
        else:
            self._radius = new_radius

    def generate_dict(self):
        """
        Generates data dict along with the Shape info
        for Cube, 
        to be used by VisualizationFrame class.
        """
        self._data_dict = {}
        self._data_dict['name'] = self._name
        self._data_dict['color'] = self._color_rgb
        self._data_dict['type'] = self.__repr__()
        self._data_dict['radius'] = self._radius
        
        return self._data_dict
