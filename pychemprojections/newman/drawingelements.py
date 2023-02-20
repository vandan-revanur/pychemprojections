import matplotlib.pyplot as plt
from dataclasses import dataclass
from math import radians
from numpy import cos, sin
from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.newman.drawingclasses import AtomInfo

logger = get_module_logger(__name__)


@dataclass
class Point:
    x: float
    y: float


def rotate_annotation(annotation_coords, angle_deg):
    rotated_annotation_coords = rotate_point_wrt_origin(annotation_coords, angle_deg)
    return rotated_annotation_coords


def rotate_point_wrt_origin(point_to_be_rotated, angle_degs):
    angle_rad = radians(angle_degs)
    point_x = point_to_be_rotated.x
    point_y = point_to_be_rotated.y
    point_x_rot = cos(angle_rad) * (point_x) - sin(angle_rad) * (point_y)
    point_y_rot = sin(angle_rad) * (point_x) + cos(angle_rad) * (point_y)

    return Point(point_x_rot, point_y_rot)


def rotate_line(point_1, point_2, angle_deg):
    point_1_rot = rotate_point_wrt_origin(point_1, angle_deg)
    point_2_rot = rotate_point_wrt_origin(point_2, angle_deg)
    xcoords = [point_1_rot.x, point_2_rot.x]
    ycoords = [point_1_rot.y, point_2_rot.y]
    return xcoords, ycoords


def draw_circle(ax, center, radius):
    bc = plt.Circle(center, radius, color="black", fill=False)
    ax.add_artist(bc)
    return ax


def add_lines(
    ax, angles_degs, pos, origin_x, origin_y, start_line_offset, end_line_offset
):
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


def add_atoms(atom_info: AtomInfo):
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
