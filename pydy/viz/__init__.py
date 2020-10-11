__all__ = []

# The following pattern is used below for importing sub-modules:
#
# 1. "from foo import *".  This imports all the names from foo.__all__ into
#    this module. But, this does not put those names into the __all__ of
#    this module. This enables "from sympy.physics.mechanics import kinematics" to
#    work.
# 2. "import foo; __all__.extend(foo.__all__)". This adds all the names in
#    foo.__all__ to the __all__ of this module. The names in __all__
#    determine which names are imported when
#    "from sympy.physics.mechanics import *" is done.

from . import visualization_frame
from .visualization_frame import *
__all__.extend(visualization_frame.__all__)

from . import shapes
from .shapes import *
__all__.extend(shapes.__all__)

from . import scene
from .scene import *
__all__.extend(scene.__all__)

from . import camera
from .camera import *
__all__.extend(camera.__all__)

from . import light
from .light import *
__all__.extend(light.__all__)

from . import server
from .server import *
__all__.extend(server.__all__)
