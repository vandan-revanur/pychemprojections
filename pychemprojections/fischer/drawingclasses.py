from dataclasses import dataclass
from typing import List


@dataclass
class DrawingInfo:
    input_smiles: str
    canvas_width_pixels: int
    canvas_height_pixels: int
    smiles_preprocessed: str
    iupac_name: str = None


@dataclass
class SingleChiralFischerNotation:
    up: str
    down: str
    left: str
    right: str


@dataclass
class MultipleChiralFischerNotation:
    up: str
    down: str
    left: List[str]
    right: List[str]
