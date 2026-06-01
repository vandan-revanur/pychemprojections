"""
Unit tests for pychemprojections.utils.smiles_manipulation_utils

This module has zero external dependencies (no RDKit, no matplotlib),
so every test here is fast and hermetic.
"""

from pychemprojections.utils.smiles_manipulation_utils import (
    get_index_of_corresponding_bracket,
)


class TestGetIndexOfCorrespondingBracket:
    """Tests for bracket-matching utility used by all projection parsers."""

    # ------------------------------------------------------------------ #
    # Happy-path: parentheses
    # ------------------------------------------------------------------ #
    def test_simple_paren_pair(self):
        assert get_index_of_corresponding_bracket("(abc)", 0) == 4

    def test_paren_at_offset(self):
        s = "XY(abc)Z"
        assert get_index_of_corresponding_bracket(s, 2) == 6

    def test_nested_parens(self):
        s = "(a(b(c)d)e)"
        assert get_index_of_corresponding_bracket(s, 0) == 10

    def test_multiple_parens_returns_first_closing(self):
        # Only the first pair: "(ab)(cd)"
        s = "(ab)(cd)"
        assert get_index_of_corresponding_bracket(s, 0) == 3

    def test_deeply_nested_parens(self):
        s = "((()))"
        assert get_index_of_corresponding_bracket(s, 0) == 5

    # ------------------------------------------------------------------ #
    # Happy-path: square brackets
    # ------------------------------------------------------------------ #
    def test_simple_square_bracket_pair(self):
        s = "[C@@H]"
        assert (
            get_index_of_corresponding_bracket(
                s, 0, start_bracket_type="[", end_bracket_type="]"
            )
            == 5
        )

    def test_square_bracket_nested(self):
        s = "[[abc]]"
        assert (
            get_index_of_corresponding_bracket(
                s, 0, start_bracket_type="[", end_bracket_type="]"
            )
            == 6
        )

    def test_square_bracket_at_offset(self):
        s = "X[C@H]Y"
        assert (
            get_index_of_corresponding_bracket(
                s, 1, start_bracket_type="[", end_bracket_type="]"
            )
            == 5
        )

    # ------------------------------------------------------------------ #
    # Real SMILES-like strings (integration-flavour)
    # ------------------------------------------------------------------ #
    def test_smiles_chiral_bracket(self):
        # "[H][C@@](F)(Cl)Br" – '[' at 3, ']' at 7
        # indices: 0=[, 1=H, 2=], 3=[, 4=C, 5=@, 6=@, 7=]
        s = "[H][C@@](F)(Cl)Br"
        assert (
            get_index_of_corresponding_bracket(
                s, 3, start_bracket_type="[", end_bracket_type="]"
            )
            == 7
        )

    def test_smiles_substituent_paren(self):
        # "(F)" starts at index 8 in "[H][C@@](F)(Cl)Br"
        # indices: ..., 7=], 8=(, 9=F, 10=)
        s = "[H][C@@](F)(Cl)Br"
        assert get_index_of_corresponding_bracket(s, 8) == 10

    # ------------------------------------------------------------------ #
    # Edge cases / error paths
    # ------------------------------------------------------------------ #
    def test_returns_minus_one_when_start_char_does_not_match(self):
        """Function contract: return -1 if char at index is not the start bracket."""
        s = "abc(def)"
        assert get_index_of_corresponding_bracket(s, 0) == -1  # 'a' != '('

    def test_returns_minus_one_for_unmatched_bracket(self):
        """Unmatched opening bracket should return -1."""
        s = "(abc"
        result = get_index_of_corresponding_bracket(s, 0)
        assert result == -1

    def test_empty_brackets(self):
        s = "()"
        assert get_index_of_corresponding_bracket(s, 0) == 1

    def test_single_char_content(self):
        s = "(X)"
        assert get_index_of_corresponding_bracket(s, 0) == 2
