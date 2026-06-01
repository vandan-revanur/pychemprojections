"""
Integration tests for pychemprojections.wedgedash.visualization

Covers plot_wedgedash_projection end-to-end with file I/O and network mocked.
The function returns the SVG string from smiles_to_svg, so we can assert on
the return value without touching the filesystem.
"""

from unittest.mock import patch

from pychemprojections.wedgedash.visualization import plot_wedgedash_projection

NO_IUPAC = patch(
    "pychemprojections.wedgedash.visualization.get_iupac_name_from_smiles",
    return_value=None,
)
NO_MAKEDIRS = patch("os.makedirs")

# Minimal fake SVG string returned by the mocked smiles_to_svg
FAKE_SVG = "<svg>test</svg>"


def _mock_svg():
    return patch(
        "pychemprojections.wedgedash.visualization.smiles_to_svg",
        return_value=FAKE_SVG,
    )


# ---------------------------------------------------------------------------
# plot_wedgedash_projection
# ---------------------------------------------------------------------------


class TestPlotWedgeDashProjection:
    CHIRAL_SMILES = "[C@@H](F)(Cl)Br"

    def test_returns_svg_string(self):
        with NO_IUPAC, NO_MAKEDIRS, _mock_svg():
            result = plot_wedgedash_projection(self.CHIRAL_SMILES)
        assert result == FAKE_SVG

    def test_runs_without_error_for_chiral_molecule(self):
        with NO_IUPAC, NO_MAKEDIRS, _mock_svg():
            plot_wedgedash_projection(self.CHIRAL_SMILES)

    def test_runs_without_error_for_two_chiral_centers(self):
        with NO_IUPAC, NO_MAKEDIRS, _mock_svg():
            plot_wedgedash_projection("[C@@H](F)(Cl)[C@H](Br)(I)")

    def test_custom_canvas_dimensions_passed_to_svg(self):
        """Canvas width/height must be forwarded to smiles_to_svg."""
        with NO_IUPAC, NO_MAKEDIRS as _md:
            with patch(
                "pychemprojections.wedgedash.visualization.smiles_to_svg",
                return_value=FAKE_SVG,
            ) as mock_svg:
                plot_wedgedash_projection(
                    self.CHIRAL_SMILES,
                    canvas_width_pixels=800,
                    canvas_height_pixels=400,
                )
                _, call_kwargs = mock_svg.call_args
                call_args = mock_svg.call_args[0]
                # smiles_to_svg(mol, path, width, height)
                assert call_args[2] == 800
                assert call_args[3] == 400

    def test_iupac_name_used_in_output_filename(self):
        """When IUPAC name is available it must appear in the output path."""
        with (
            patch(
                "pychemprojections.wedgedash.visualization.get_iupac_name_from_smiles",
                return_value="bromochlorofluoromethane",
            ),
            NO_MAKEDIRS,
        ):
            with patch(
                "pychemprojections.wedgedash.visualization.smiles_to_svg",
                return_value=FAKE_SVG,
            ) as mock_svg:
                plot_wedgedash_projection(self.CHIRAL_SMILES)
                output_path = mock_svg.call_args[0][1]
                assert "bromochlorofluoromethane" in output_path

    def test_fallback_filename_when_no_iupac(self):
        """Without an IUPAC name a default filename must be used."""
        with NO_IUPAC, NO_MAKEDIRS:
            with patch(
                "pychemprojections.wedgedash.visualization.smiles_to_svg",
                return_value=FAKE_SVG,
            ) as mock_svg:
                plot_wedgedash_projection(self.CHIRAL_SMILES)
                output_path = mock_svg.call_args[0][1]
                assert "output_wedgedash" in output_path
