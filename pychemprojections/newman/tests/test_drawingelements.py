"""
Unit tests for pychemprojections.newman.drawingelements

Geometry functions (rotate_point_wrt_origin, rotate_line, rotate_annotation)
are pure math – no external dependencies.  Drawing functions (draw_circle,
add_lines, add_atoms) require a matplotlib Axes object and are tested via
integration smoke tests that verify the Axes is returned without errors.
"""

import math
import matplotlib

matplotlib.use("Agg")  # non-interactive backend – safe for CI
import matplotlib.pyplot as plt

from pychemprojections.newman.drawingelements import (
    Point,
    rotate_point_wrt_origin,
    rotate_annotation,
    rotate_line,
    draw_circle,
    add_lines,
    add_atoms,
)
from pychemprojections.newman.drawingclasses import AtomInfo


TOLERANCE = 1e-9


# ---------------------------------------------------------------------------
# rotate_point_wrt_origin
# ---------------------------------------------------------------------------
class TestRotatePointWrtOrigin:
    def test_rotate_0_degrees_is_identity(self):
        p = Point(3.0, 4.0)
        result = rotate_point_wrt_origin(p, 0)
        assert abs(result.x - 3.0) < TOLERANCE
        assert abs(result.y - 4.0) < TOLERANCE

    def test_rotate_90_degrees(self):
        # (1,0) → (0,1) under 90° CCW
        p = Point(1.0, 0.0)
        result = rotate_point_wrt_origin(p, 90)
        assert abs(result.x - 0.0) < TOLERANCE
        assert abs(result.y - 1.0) < TOLERANCE

    def test_rotate_180_degrees(self):
        # (1,0) → (-1,0) under 180°
        p = Point(1.0, 0.0)
        result = rotate_point_wrt_origin(p, 180)
        assert abs(result.x - (-1.0)) < TOLERANCE
        assert abs(result.y - 0.0) < TOLERANCE

    def test_rotate_270_degrees(self):
        # (0,1) → (1,0) under 270° CCW  (= -90°)
        p = Point(0.0, 1.0)
        result = rotate_point_wrt_origin(p, 270)
        assert abs(result.x - 1.0) < TOLERANCE
        assert abs(result.y - 0.0) < TOLERANCE

    def test_rotate_360_degrees_is_identity(self):
        p = Point(2.5, -3.1)
        result = rotate_point_wrt_origin(p, 360)
        assert abs(result.x - p.x) < TOLERANCE
        assert abs(result.y - p.y) < TOLERANCE

    def test_returns_point_instance(self):
        result = rotate_point_wrt_origin(Point(1, 0), 45)
        assert isinstance(result, Point)

    def test_rotate_45_degrees(self):
        p = Point(1.0, 0.0)
        result = rotate_point_wrt_origin(p, 45)
        expected = math.sqrt(2) / 2
        assert abs(result.x - expected) < TOLERANCE
        assert abs(result.y - expected) < TOLERANCE


# ---------------------------------------------------------------------------
# rotate_annotation (delegates to rotate_point_wrt_origin)
# ---------------------------------------------------------------------------
class TestRotateAnnotation:
    def test_same_result_as_rotate_point(self):
        p = Point(1.0, 0.0)
        angle = 60.0
        from_annotation = rotate_annotation(p, angle)
        from_point = rotate_point_wrt_origin(p, angle)
        assert abs(from_annotation.x - from_point.x) < TOLERANCE
        assert abs(from_annotation.y - from_point.y) < TOLERANCE


# ---------------------------------------------------------------------------
# rotate_line
# ---------------------------------------------------------------------------
class TestRotateLine:
    def test_returns_tuple_of_two_lists(self):
        xc, yc = rotate_line(Point(0, 0), Point(0, 1), 0)
        assert isinstance(xc, list) and isinstance(yc, list)

    def test_each_list_has_two_elements(self):
        xc, yc = rotate_line(Point(0, 0), Point(1, 0), 90)
        assert len(xc) == 2 and len(yc) == 2

    def test_zero_rotation_preserves_coords(self):
        xc, yc = rotate_line(Point(1, 2), Point(3, 4), 0)
        assert abs(xc[0] - 1) < TOLERANCE
        assert abs(yc[0] - 2) < TOLERANCE
        assert abs(xc[1] - 3) < TOLERANCE
        assert abs(yc[1] - 4) < TOLERANCE

    def test_90_degree_rotation_vertical_line(self):
        # Vertical line (0,0)→(0,1) rotated 90° → horizontal (0,0)→(-1,0)
        xc, yc = rotate_line(Point(0, 0), Point(0, 1), 90)
        assert abs(xc[0] - 0.0) < TOLERANCE
        assert abs(yc[0] - 0.0) < TOLERANCE
        assert abs(xc[1] - (-1.0)) < TOLERANCE
        assert abs(yc[1] - 0.0) < TOLERANCE


# ---------------------------------------------------------------------------
# draw_circle  (smoke test – just verify Axes returned without error)
# ---------------------------------------------------------------------------
class TestDrawCircle:
    def test_returns_axes(self):
        fig, ax = plt.subplots()
        result = draw_circle(ax, (0.5, 0.5), 0.3)
        assert result is ax
        plt.close(fig)

    def test_circle_added_to_axes_artists(self):
        fig, ax = plt.subplots()
        draw_circle(ax, (0, 0), 0.1)
        # At least one artist (the circle) should have been added
        assert len(ax.get_children()) > 0
        plt.close(fig)


# ---------------------------------------------------------------------------
# add_lines  (smoke test)
# ---------------------------------------------------------------------------
class TestAddLines:
    def _make_atom_info(self, ax, pos="front"):
        return AtomInfo(
            ax=ax,
            atoms=["CH3", "Br", "Cl"],
            angles=[0, 120, 240],
            fontsize=10,
            origin_x=0.5,
            origin_y=0.5,
            annotation_offset=0.2,
            pos=pos,
        )

    def test_front_returns_axes(self):
        fig, ax = plt.subplots()
        atom_info = self._make_atom_info(ax, pos="front")
        result = add_lines(atom_info, start_line_offset=0.05, end_line_offset=0.3)
        assert result is ax
        plt.close(fig)

    def test_rear_returns_axes(self):
        fig, ax = plt.subplots()
        atom_info = self._make_atom_info(ax, pos="rear")
        result = add_lines(atom_info, start_line_offset=0.05, end_line_offset=0.3)
        assert result is ax
        plt.close(fig)


# ---------------------------------------------------------------------------
# add_atoms  (smoke test)
# ---------------------------------------------------------------------------
class TestAddAtoms:
    def test_returns_axes(self):
        fig, ax = plt.subplots()
        atom_info = AtomInfo(
            ax=ax,
            atoms=["$CH_{3}$", "$Br$", "$Cl$"],
            angles=[0, 120, 240],
            fontsize=10,
            origin_x=0.5,
            origin_y=0.5,
            annotation_offset=0.2,
            pos="front",
        )
        result = add_atoms(atom_info)
        assert result is ax
        plt.close(fig)
