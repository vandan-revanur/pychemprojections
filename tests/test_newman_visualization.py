"""
Integration tests for pychemprojections.newman.visualization

Covers:
  - get_newman_drawing_info  (full pipeline, network mocked)
  - get_front_and_back_groups_for_plotting  (pure logic)
  - plot_newman_projection  (smoke – file I/O mocked)
"""

import matplotlib

matplotlib.use("Agg")

from unittest.mock import patch
import matplotlib.pyplot as plt

from pychemprojections.newman.visualization import (
    get_newman_drawing_info,
    get_front_and_back_groups_for_plotting,
    plot_newman_projection,
)
from pychemprojections.newman.drawingclasses import NewmanDrawingInfo

NO_IUPAC = patch(
    "pychemprojections.newman.visualization.get_iupac_name_from_smiles",
    return_value=None,
)
NO_SAVE = patch("matplotlib.figure.Figure.savefig")
NO_MAKEDIRS = patch("os.makedirs")
NO_SHOW = patch("matplotlib.pyplot.show")


# ---------------------------------------------------------------------------
# get_front_and_back_groups_for_plotting  (pure)
# ---------------------------------------------------------------------------


class TestGetFrontAndBackGroupsForPlotting:
    # Inputs are raw condensed strings; the function wraps them in $...$
    # via prepare_string_for_newman_projection_plot.
    GROUPS = ["CH3", "F", "Cl", "Br", "I", "H"]

    def test_first_three_are_front(self):
        front, _ = get_front_and_back_groups_for_plotting(self.GROUPS)
        assert len(front) == 3

    def test_last_three_are_rear(self):
        _, rear = get_front_and_back_groups_for_plotting(self.GROUPS)
        assert len(rear) == 3

    def test_returns_two_lists_of_three(self):
        front, rear = get_front_and_back_groups_for_plotting(self.GROUPS)
        assert len(front) == 3
        assert len(rear) == 3

    def test_output_strings_wrapped_in_dollars(self):
        front, rear = get_front_and_back_groups_for_plotting(self.GROUPS)
        for label in front + rear:
            assert label.startswith("$") and label.endswith("$")

    def test_order_is_preserved(self):
        """First input becomes first front atom, fourth becomes first rear atom."""
        front, rear = get_front_and_back_groups_for_plotting(self.GROUPS)
        # CH3 → $CH_{3}$, Br → $Br$
        assert "CH" in front[0]
        assert "Br" in rear[0]


# ---------------------------------------------------------------------------
# get_newman_drawing_info  (full pipeline)
# ---------------------------------------------------------------------------


class TestGetNewmanDrawingInfo:
    # 1-bromobutane: clear C-C bond to examine
    SMILES = "CCCCBr"

    def test_returns_newman_drawing_info(self):
        with NO_IUPAC:
            result = get_newman_drawing_info(self.SMILES)
        assert isinstance(result, NewmanDrawingInfo)

    def test_front_atoms_has_three_entries(self):
        with NO_IUPAC:
            result = get_newman_drawing_info(self.SMILES)
        assert len(result.front_atoms) == 3

    def test_rear_atoms_has_three_entries(self):
        with NO_IUPAC:
            result = get_newman_drawing_info(self.SMILES)
        assert len(result.rear_atoms) == 3

    def test_all_atom_labels_are_dollar_wrapped(self):
        with NO_IUPAC:
            result = get_newman_drawing_info(self.SMILES)
        for label in result.front_atoms + result.rear_atoms:
            assert label.startswith("$") and label.endswith("$"), (
                f"Label '{label}' is not wrapped in $ signs"
            )

    def test_explicit_bond_ids_respected(self):
        """When carbon_ids_bond_to_examine is given, the function must not crash."""
        with NO_IUPAC:
            result = get_newman_drawing_info(
                self.SMILES, carbon_ids_bond_to_examine=(0, 1)
            )
        assert isinstance(result, NewmanDrawingInfo)

    def test_canvas_dimensions_stored(self):
        with NO_IUPAC:
            result = get_newman_drawing_info(
                self.SMILES, canvas_width_pixels=800, canvas_height_pixels=400
            )
        assert result.canvas_width_pixels == 800
        assert result.canvas_height_pixels == 400

    def test_iupac_name_stored_when_provided(self):
        with patch(
            "pychemprojections.newman.visualization.get_iupac_name_from_smiles",
            return_value="1-bromobutane",
        ):
            result = get_newman_drawing_info(self.SMILES)
        assert result.iupac_name == "1-bromobutane"


# ---------------------------------------------------------------------------
# plot_newman_projection  (smoke)
# ---------------------------------------------------------------------------


class TestPlotNewmanProjection:
    def _make_info(self, iupac_name=None):
        return NewmanDrawingInfo(
            front_atoms=["$CH_{3}$", "$F$", "$Cl$"],
            rear_atoms=["$Br$", "$H$", "$H$"],
            canvas_width_pixels=600,
            canvas_height_pixels=600,
            iupac_name=iupac_name,
        )

    def test_runs_without_error_no_iupac(self):
        with NO_SAVE, NO_MAKEDIRS, NO_SHOW:
            plot_newman_projection(self._make_info())
        plt.close("all")

    def test_runs_without_error_with_iupac(self):
        with NO_SAVE, NO_MAKEDIRS, NO_SHOW:
            plot_newman_projection(self._make_info(iupac_name="1-bromobutane"))
        plt.close("all")
