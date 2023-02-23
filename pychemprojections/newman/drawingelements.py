import matplotlib.pyplot as plt
from dataclasses import dataclass
from math import radians
from numpy import cos, sin
from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.newman.drawingclasses import AtomInfo
from matplotlib.axes import Axes
from typing import List, Tuple

logger = get_module_logger(__name__)


@dataclass
class Point:
    x: float
    y: float


def rotate_annotation(annotation_coords: Point, angle_deg: float):
    """
    Rotate the coordinates of the annotation labels of the atoms

    Parameters
    ----------
    annotation_coords :Point
        Point containing coordinates of the annotation labels of the atom

    angle_deg : float
        Angle in degrees to rotate the point

    Returns
    -------
    Point
    Rotated Annotation label point coordinates

    """
    rotated_annotation_coords = rotate_point_wrt_origin(annotation_coords, angle_deg)
    return rotated_annotation_coords


def rotate_point_wrt_origin(point_to_be_rotated: Point, angle_degs: float) -> Point:
    """
    Rotate a point with respect to the origin

    Parameters
    ----------
    point_to_be_rotated : Point
        Point to be rotated

    angle_degs : float
        Angle in degrees to rotate the point

    Returns
    -------
    Point
    Point rotated about the origin

    """
    angle_rad = radians(angle_degs)
    point_x = point_to_be_rotated.x
    point_y = point_to_be_rotated.y
    point_x_rot = cos(angle_rad) * (point_x) - sin(angle_rad) * (point_y)
    point_y_rot = sin(angle_rad) * (point_x) + cos(angle_rad) * (point_y)

    return Point(point_x_rot, point_y_rot)


def rotate_line(point_1: Point, point_2, angle_deg) -> Tuple[List[float], List[float]]:
    """
    Rotate a line

    Parameters
    ----------
    point_1 : Point
        Start point of the line

    point_2 : Point
        End point of the line

    angle_deg : float
        Angle to rotate the line

    Returns
    -------
    Tuple[List[float], List[float]]
    A tuple of the x coordinates and y coordinates after rotation

    """
    point_1_rot = rotate_point_wrt_origin(point_1, angle_deg)
    point_2_rot = rotate_point_wrt_origin(point_2, angle_deg)
    xcoords = [point_1_rot.x, point_2_rot.x]
    ycoords = [point_1_rot.y, point_2_rot.y]
    return xcoords, ycoords


def draw_circle(ax: Axes, center: Tuple[float, float], radius: float) -> Axes:
    """
    Draw the circle in the Newman projection

    Parameters
    ----------
    ax : Axes
        A matplotlib Axes object

    center : Tuple[float,float]
        Coordinates of the center of the circle

    radius : float
        Radius of the circle

    Returns
    -------
    Axes
    A matplotlib Axes object

    """
    bc = plt.Circle(center, radius, color="black", fill=False)
    ax.add_artist(bc)
    return ax


def add_lines(
    atom_info: AtomInfo, start_line_offset: float, end_line_offset: float
) -> Axes:
    """
    Add lines to the newman projection diagram

    Parameters
    ----------
    atom_info : AtomInfo
        Atom related information such as front and back atoms, front and back angles, axes, fontsize, origin coordinates,
        and position whether it is front or rear atoms

    start_line_offset : float
        Coordinate offset to add at the point where the lines are supposed to start

    end_line_offset : float
        Coordinate offset to add at the point where the lines are supposed to end

    Returns
    -------
    Axes
    A matplotlib Axes object

    """
    origin_x = atom_info.origin_x
    origin_y = atom_info.origin_y
    pos = atom_info.pos
    angles_degs = atom_info.angles
    ax = atom_info.ax

    line_end_coords = Point(origin_x, origin_y + end_line_offset)
    if pos == "rear":
        line_start_coords = Point(origin_x, origin_y + start_line_offset)
    else:
        line_start_coords = Point(origin_x, origin_y)

    for line_idx, angle in enumerate(angles_degs):
        logger.debug(f"rotating line  {line_idx} in {pos} by {angle} degrees")
        xcoords, ycoords = rotate_line(line_start_coords, line_end_coords, angle)
        ax.plot(xcoords, ycoords, color="black")
    return ax


def add_atoms(atom_info: AtomInfo) -> Axes:
    """
    Add atoms to the newman projection diagram

    Parameters
    ----------
    atom_info : AtomInfo
        Atom related information such as front and back atoms, front and back angles etc

    Returns
    -------
    Axes
    A matplotlib Axes object

    """
    ax = atom_info.ax
    labels = atom_info.atoms
    angles_degs = atom_info.angles
    fontsize = atom_info.fontsize
    origin_x = atom_info.origin_x
    origin_y = atom_info.origin_y
    annotation_offset = atom_info.annotation_offset
    pos = atom_info.pos

    annotation_start_coords = Point(origin_x, origin_y + annotation_offset)

    logger.info(f"atom labels in the {pos}: {labels}")
    for atom_idx, angle in enumerate(angles_degs):
        atom_label = labels[atom_idx]
        logger.debug(
            f"rotating atom {atom_idx} by {angle} degrees and assigning label: {atom_label}"
        )
        annotation_coords = rotate_annotation(annotation_start_coords, angle)

        ax.annotate(
            atom_label,
            [annotation_coords.x, annotation_coords.y],
            ha="center",
            va="center",
            color="black",
            fontweight="bold",
            fontsize=fontsize,
        )
    return ax
