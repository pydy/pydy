#importing matplotlib's colorconverter here,
#module can be copied later in our source.
__all__ = ['Shape']


from matplotlib.colors import ColorConverter


convert = ColorConverter()

class Shape(object):
    """
    A Shape. It is a superclass for more general shapes like
    Mesh, Cylinder, Sphere etc. 
    
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
    >>> s.name()
    'UnNamed'
    >>> s.color()
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color()
    'blue'
    >>> a = Shape('my-shape2', 'red')
    >>> a.name()
    'my-shape2'
    >>> a.color()
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

    def _generate_dict(self):
        """
        Generates data dict along with the Shape info
        to be used by VisualizationFrame class.
        """
        
    def __str__(self):
        return 'Shape ' + self._name + ' color:' + self._color
    
    def __repr__(self):
        return 'Shape'    

    def generate_visualization_frame(self, reference_frame, origin, \
                                                      name='UnNamed'):
        ## I am not sure whether we want to include this method.
        ## Or we want to generate VisualizationFrames from Scene class 
        #Instead. 
        
        """
        creates a VisualizationFrame with origin as the center
        of this shape.
        
        Parameters
        ==========
        reference_frame: an Instance of mechanics.ReferenceFrame
        The frame of reference the VisualizationFrame is oriented in.
        origin: an Instance of mechanics.Point
        The origin point of the VisualizationFrame
        
        Examples
        ========

        >>> from pydy_viz import Shape
        >>> from sympy.physics.mechanics import ReferenceFrame, Point
        >>> s = Shape()
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> f = s.generate_visualization_frame(I,O,name='frame1')
        
        """
        ##Writing the general concept for initializing a visualization
        #frame. It might need to change with the addition of 
        #VisualizationFrame class
        
        if not isinstance(name, str):
            raise TypeError('name should be a valid str object')
        else:
            _frame = VisualizationFrame(name)
        
            
        if not isinstance(reference_frame, ReferenceFrame):
            raise TypeError('%s should be a valid \
                                     ReferenceFrame'%reference_frame)        
        else:
            _frame.reference_frame = reference_frame                             

        if not isinstance(origin, Point):
            raise TypeError('%s should be a valid \
                                     Point'%origin)        
        else:
            _frame.origin = origin       
        
        return _frame
            
class Cube(object):
    """
    A Cube. This class generates a Cube, with given length of side, 
    and color.Default color is Grey.
    
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
    >>> s.name()
    'UnNamed'
    >>> s.color()
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>>s.length()
    10
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length()
    12
    >>> a = Cube(10, 'my-shape2', 'red')
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
                                            'length: ' + self._length
    
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

class Cylinder(object):
    """
    A Cylinder. This class generates a Cylinder with given length,
    radius, and color. Default color is Grey.
    
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
    >>> s = Cylinder(10, 5)
    >>> s.name()
    'UnNamed'
    >>> s.color()
    'grey'
    >>> s.color_in_rgb()
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>> s.length()
    10
    >>> s.radius()
    5
    >>>#These can be changed later too ..
    >>> s.name = 'my-shape1'
    >>> s.name
    'my-shape1'
    >>> s.color = 'blue'
    >>> s.color
    'blue'
    >>> s.length = 12
    >>> s.length()
    12
    >>> s.radius = 6
    >>> s.radius()
    6
    >>> a = Cylinder(10, 5, 'my-shape2', 'red')
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
                                        'length: ' + self._length + \
                                        'radius: ' + self._radius
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
