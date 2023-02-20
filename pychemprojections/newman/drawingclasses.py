from dataclasses import dataclass
from typing import List
from matplotlib.axes import Axes


@dataclass()
class NewmanDrawingInfo:
    rear_atoms: List[str]
    front_atoms: List[str]
    canvas_width_pixels: int
    canvas_height_pixels: int
    iupac_name: str


@dataclass()
class AtomInfo:
    ax: Axes
    atoms: List[str]
    angles: List[int]
    fontsize: int
    origin_x: float
    origin_y: float
    annotation_offset: float
    pos: str
