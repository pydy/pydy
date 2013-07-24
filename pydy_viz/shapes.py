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

    >>> from pydy_viz import Shape
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
    >>> a = Shape('my-shape2','red')
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

    def _generate_dict(self):
        """
        Generates data dict along with the Shape info
        to be used by VisualizationFrame class.
        """
        
    def __str__(self):
        return 'Shape ' + self._name + ' color:' + self._color
    
    def __repr__(self):
        return 'Shape'    

    def generate_visualization_frame(self, reference_frame, origin,\
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
        
        if not isinstance(name,str):
            raise TypeError('name should be a valid str object')
        else:
            _frame = VisualizationFrame(name)
        
            
        if not isinstance(reference_frame,ReferenceFrame):
            raise TypeError('%s should be a valid \
                                     ReferenceFrame')%(reference_frame)        
        else:
            _frame.reference_frame = reference_frame                             

        if not isinstance(origin,Point):
            raise TypeError('%s should be a valid \
                                     Point')%(point)        
        else:
            _frame.origin = origin       
        
        return _frame
            
