"""
Unit tests for pychemprojections.fischer.stringmanipulation

get_width_of_text_string and prepare_strings_for_fischer_projection_plot both
spin up a matplotlib figure; they are tested with a small canvas so tests run
quickly without producing files.

split_into_chunks and get_condensed_smiles_of_all_substituents are pure
logic and are tested exhaustively.
"""

from pychemprojections.fischer.stringmanipulation import (
    split_into_chunks,
    get_condensed_smiles_of_all_substituents,
    prepare_strings_for_fischer_projection_plot,
)


# ---------------------------------------------------------------------------
# split_into_chunks
# ---------------------------------------------------------------------------
class TestSplitIntoChunks:
    """
    split_into_chunks(inp_str, text_width_pixels_coords, threshold_pixels)
    Wraps a string whose total width exceeds threshold by inserting "-\n-$" delimiters.
    """

    def test_short_string_not_split(self):
        # text_width < threshold → no split, input returned directly
        result = split_into_chunks("CH_3_", 50, 200)
        # No line-break separator added
        assert "-\n-$" not in result

    def test_split_result_contains_separator(self):
        # split_into_chunks expects '_'-marked strings produced upstream.
        # Ratio: text_width=120, len=12 → pixels_per_char=10,
        #        threshold=40 → threshold_chars=4 → splits occur.
        inp = "C_2_H_5_O_H_"
        result = split_into_chunks(inp, 120, 40)
        assert "-\n-$" in result

    def test_trailing_separator_stripped(self):
        """The function removes a trailing '-\n-$' if present."""
        # We can't trivially produce the exact edge case without exact pixel
        # arithmetic, so verify the invariant on any result.
        inp = "C_2_H_5_"
        result = split_into_chunks(inp, 100, 20)
        assert not result.endswith("-\n-$")

    def test_returns_string(self):
        assert isinstance(split_into_chunks("CH_3_", 100, 200), str)


# ---------------------------------------------------------------------------
# get_condensed_smiles_of_all_substituents
# ---------------------------------------------------------------------------
class TestGetCondensedSmilesOfAllSubstituents:
    """
    Delegates to smiles_to_condensed_form(cleanup_Hs(s)) for each substituent.
    Test end-to-end condensed output.
    """

    def test_returns_list_of_same_length(self):
        subs = ["[H]CC", "Br", "Cl"]
        result = get_condensed_smiles_of_all_substituents(subs)
        assert len(result) == len(subs)

    def test_explicit_H_brackets_removed(self):
        result = get_condensed_smiles_of_all_substituents(["[H]CC"])
        assert "[H]" not in result[0]

    def test_plain_substituent_passes_through(self):
        result = get_condensed_smiles_of_all_substituents(["Br"])
        assert result[0] == "Br"

    def test_all_results_are_strings(self):
        subs = ["[H]CC", "F", "Br", "Cl"]
        result = get_condensed_smiles_of_all_substituents(subs)
        assert all(isinstance(s, str) for s in result)

    def test_empty_list_returns_empty_list(self):
        assert get_condensed_smiles_of_all_substituents([]) == []


# ---------------------------------------------------------------------------
# prepare_strings_for_fischer_projection_plot  (integration, small canvas)
# ---------------------------------------------------------------------------
class TestPrepareStringsForFischerProjectionPlot:
    CANVAS_W = 600
    CANVAS_H = 400

    def test_result_wrapped_in_dollar_signs(self):
        result = prepare_strings_for_fischer_projection_plot(
            "CH3", self.CANVAS_W, self.CANVAS_H
        )
        assert result.startswith("$") and result.endswith("$")

    def test_numbers_preceded_by_underscore(self):
        result = prepare_strings_for_fischer_projection_plot(
            "C2H5", self.CANVAS_W, self.CANVAS_H
        )
        # Numbers must have an underscore prefix inserted by the function
        assert "_2" in result or "_5" in result

    def test_returns_string(self):
        result = prepare_strings_for_fischer_projection_plot(
            "Br", self.CANVAS_W, self.CANVAS_H
        )
        assert isinstance(result, str)

    def test_atom_without_numbers_no_underscore_added(self):
        result = prepare_strings_for_fischer_projection_plot(
            "Br", self.CANVAS_W, self.CANVAS_H
        )
        # "Br" has no digits; no underscores should be inserted for subscripts
        assert "_" not in result
