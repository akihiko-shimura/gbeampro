__version__ = "2.2.0"

from .beam import GaussBeam
from .elements import (
    Element,
    Propagation,
    ThinLens,
    Interface,
    InterfaceCurved,
    CurvedMirrorTan,
    CurvedMirrorSag,
)
from .system import OpticalSystem
from . import analysis, plot, optimize
