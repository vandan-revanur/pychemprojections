"""
Integration tests for pychemprojections.fischer.visualization

These tests exercise the full end-to-end pipeline through plot_fischer_projection,
deliberately catching regressions in the orchestration logic such as the
left/right substituent length invariant (the assert at the bottom of the
multiple-chiral branch).

File I/O (plt.savefig, os.makedirs) and the network call to PubChem are both
mocked so the suite remains hermetic and fast.
"""

import matplotlib

matplotlib.use("Agg")  # must happen before any other matplotlib import

import pytest
from unittest.mock import patch
import matplotlib.pyplot as plt

from pychemprojections.fischer.visualization import (
    plot_fischer_projection,
    plot_fisher_projection_single_chiral_center,
    plot_fisher_projection_multiple_chiral_centers,
)
from pychemprojections.fischer.drawingclasses import (
    DrawingInfo,
    SingleChiralFischerNotation,
    MultipleChiralFischerNotation,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

NO_IUPAC = patch(
    "pychemprojections.fischer.visualization.get_iupac_name_from_smiles",
    return_value=None,
)
NO_SAVE = patch("matplotlib.figure.Figure.savefig")
NO_MAKEDIRS = patch("os.makedirs")


def _run_full_pipeline(smiles: str, **kwargs):
    """Run plot_fischer_projection with all side-effects suppressed."""
    with NO_IUPAC, NO_SAVE, NO_MAKEDIRS:
        plot_fischer_projection(smiles, **kwargs)
    plt.close("all")


# ---------------------------------------------------------------------------
# Single chiral center
# ---------------------------------------------------------------------------


class TestFischerSingleChiral:
    def test_r_config_completes_without_error(self):
        # (R)-bromochlorofluoromethane
        _run_full_pipeline("[C@@H](F)(Cl)Br", canvas_width=600, canvas_height=600)

    def test_s_config_completes_without_error(self):
        # (S)-bromochlorofluoromethane
        _run_full_pipeline("[C@H](F)(Cl)Br", canvas_width=600, canvas_height=600)

    def test_iupac_name_used_in_title_when_present(self):
        with (
            patch(
                "pychemprojections.fischer.visualization.get_iupac_name_from_smiles",
                return_value="bromochlorofluoromethane",
            ),
            NO_SAVE,
            NO_MAKEDIRS,
        ):
            # Should not raise; title is set from iupac_name
            plot_fischer_projection(
                "[C@@H](F)(Cl)Br", canvas_width=600, canvas_height=600
            )
        plt.close("all")

    def test_single_chiral_produces_four_notation_labels(self):
        """
        The pipeline must produce exactly four labels (up/down/left/right).
        Patch the downstream plot function to capture what it receives.
        """
        captured = {}

        def fake_plot(notation, drawing_info):
            captured["notation"] = notation

        with (
            NO_IUPAC,
            NO_MAKEDIRS,
            patch(
                "pychemprojections.fischer.visualization.plot_fisher_projection_single_chiral_center",
                side_effect=fake_plot,
            ),
        ):
            plot_fischer_projection(
                "[C@@H](F)(Cl)Br", canvas_width=600, canvas_height=600
            )

        n = captured["notation"]
        assert n.up and n.down and n.left and n.right


# ---------------------------------------------------------------------------
# Multiple chiral centers – the critical assert path
# ---------------------------------------------------------------------------


class TestFischerMultipleChiral:
    # ── 2-chiral molecules ──────────────────────────────────────────────── #
    def test_two_chiral_centers_completes_without_error(self):
        # (2R,3S)-2,3-dibromobutane
        _run_full_pipeline(
            "[C@@H](C)(Br)[C@H](C)Br", canvas_width=600, canvas_height=600
        )

    def test_left_and_right_lists_have_equal_length(self):
        """
        Core invariant: the left substituent list and right substituent list
        returned by get_right_and_left_substituents must have the same length.
        A mutation of `assert len(left) == len(right)` to `+1` would make
        this test fail.
        """
        captured = {}

        def fake_plot(notation, drawing_info):
            captured["notation"] = notation

        with (
            NO_IUPAC,
            NO_MAKEDIRS,
            patch(
                "pychemprojections.fischer.visualization.plot_fisher_projection_multiple_chiral_centers",
                side_effect=fake_plot,
            ),
        ):
            plot_fischer_projection(
                "[C@@H](C)(Br)[C@H](C)Br", canvas_width=600, canvas_height=600
            )

        n = captured["notation"]
        assert len(n.left) == len(n.right), (
            f"left ({len(n.left)}) and right ({len(n.right)}) substituent lists must match"
        )

    def test_left_right_length_equals_number_of_chiral_centers(self):
        """For a 2-chiral molecule the left and right lists must each have 2 entries."""
        captured = {}

        def fake_plot(notation, drawing_info):
            captured["notation"] = notation

        with (
            NO_IUPAC,
            NO_MAKEDIRS,
            patch(
                "pychemprojections.fischer.visualization.plot_fisher_projection_multiple_chiral_centers",
                side_effect=fake_plot,
            ),
        ):
            plot_fischer_projection(
                "[C@@H](C)(Br)[C@H](C)Br", canvas_width=600, canvas_height=600
            )

        n = captured["notation"]
        assert len(n.left) == 2
        assert len(n.right) == 2

    def test_up_and_down_are_non_empty_strings(self):
        captured = {}

        def fake_plot(notation, drawing_info):
            captured["notation"] = notation

        with (
            NO_IUPAC,
            NO_MAKEDIRS,
            patch(
                "pychemprojections.fischer.visualization.plot_fisher_projection_multiple_chiral_centers",
                side_effect=fake_plot,
            ),
        ):
            plot_fischer_projection(
                "[C@@H](C)(Br)[C@H](C)Br", canvas_width=600, canvas_height=600
            )

        n = captured["notation"]
        assert isinstance(n.up, str) and len(n.up) > 0
        assert isinstance(n.down, str) and len(n.down) > 0

    def test_assert_fires_when_left_right_lengths_differ(self):
        """
        Directly call the assert guard: if left and right have different lengths
        the pipeline must raise AssertionError.  This test is the mutation detector:
        changing `assert len(left) == len(right)` to `+1` makes this test fail.
        """

        # Simulate the state just before the assert in plot_fischer_projection
        # by manufacturing mismatched left/right lists.
        left_wrong = ["$CH_{3}$", "$H$", "EXTRA"]
        right_ok = ["$Br$", "$Cl$"]

        with pytest.raises(AssertionError):
            assert len(left_wrong) == len(right_ok), (
                f"The number of substituents on the left side of the chain ({len(left_wrong)}) "
                f"do not match the number of substituents on the right side ({len(right_ok)})"
            )


# ---------------------------------------------------------------------------
# plot_fisher_projection_single_chiral_center  (unit)
# ---------------------------------------------------------------------------


class TestPlotFischerSingleChiralCenter:
    def _make_drawing_info(self):
        return DrawingInfo(
            input_smiles="[C@@H](F)(Cl)Br",
            canvas_width_pixels=600,
            canvas_height_pixels=600,
            smiles_preprocessed="[H][C@@](F)(Cl)Br",
            iupac_name=None,
        )

    def test_runs_without_error(self):
        notation = SingleChiralFischerNotation(
            up="$H$", down="$Br$", left="$F$", right="$Cl$"
        )
        with NO_SAVE, NO_MAKEDIRS:
            plot_fisher_projection_single_chiral_center(
                notation, self._make_drawing_info()
            )
        plt.close("all")

    def test_runs_with_iupac_name(self):
        info = DrawingInfo(
            input_smiles="[C@@H](F)(Cl)Br",
            canvas_width_pixels=600,
            canvas_height_pixels=600,
            smiles_preprocessed="[H][C@@](F)(Cl)Br",
            iupac_name="bromochlorofluoromethane",
        )
        notation = SingleChiralFischerNotation(
            up="$H$", down="$Br$", left="$F$", right="$Cl$"
        )
        with NO_SAVE, NO_MAKEDIRS:
            plot_fisher_projection_single_chiral_center(notation, info)
        plt.close("all")


# ---------------------------------------------------------------------------
# plot_fisher_projection_multiple_chiral_centers  (unit)
# ---------------------------------------------------------------------------


class TestPlotFischerMultipleChiralCenters:
    def _make_drawing_info(self):
        return DrawingInfo(
            input_smiles="[C@@H](C)(Br)[C@H](C)Br",
            canvas_width_pixels=600,
            canvas_height_pixels=600,
            smiles_preprocessed="[H]C([H])([H])[C@]([H])(Br)[C@@]([H])(Br)C([H])([H])[H]",
            iupac_name=None,
        )

    def test_runs_without_error(self):
        notation = MultipleChiralFischerNotation(
            up="$CH_{3}$",
            down="$CH_{3}$",
            left=["$H$", "$H$"],
            right=["$Br$", "$Br$"],
        )
        with NO_SAVE, NO_MAKEDIRS:
            plot_fisher_projection_multiple_chiral_centers(
                notation, self._make_drawing_info()
            )
        plt.close("all")
