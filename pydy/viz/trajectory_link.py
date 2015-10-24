from IPython.html.widgets import Widget
from IPython.utils.traitlets import CFloat, Instance, List, Tuple, Unicode
from pythreejs import Mesh


def vector3(trait_type=CFloat, default=None, **kwargs):
    if default is None:
        default=[0, 0, 0]
    return List(trait_type, default_value=default,
                minlen=3, maxlen=3, allow_none=False, **kwargs)


def vector4(trait_type=CFloat, default=None, **kwargs):
    if default is None:
        default=[0, 0, 0, 0]
    return List(trait_type, default_value=default,
                minlen=4, maxlen=4, allow_none=False, **kwargs)


class TrajectoryLink(Widget):
    _model_module = Unicode('nbextensions/pydyviz/pydyviz', sync=True)
    _model_name = Unicode('TrajectoryLinkModel', sync=True)

    source = Instance(Widget, sync=True)
    target = Instance(Mesh, sync=True)
    position = List(vector3(), default_value=[[0, 0, 0]], sync=True)
    quaternion = List(vector4(), default_value=[[0, 0, 0, 0]], sync=True)

    """TrajectoryLink Widget

    source: a bounded IntWidget or bounded FloatWidget
    target: a Mesh with position and quaternion trajectories
    """
    # for compatibility with traitlet links
    def unlink(self):
        self.close()


def trajectory_link(int_widget, mesh):
    return TrajectoryLink(source=int_widget, target=mesh)
