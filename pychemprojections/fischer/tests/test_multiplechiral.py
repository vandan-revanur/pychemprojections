"""
Unit tests for pychemprojections.fischer.multiplechiral

Fixture molecule: (2R,3S)-2,3-dibromobutane
Prepared SMILES: "[H]C([H])([H])[C@]([H])(Br)[C@@]([H])(Br)C([H])([H])[H]"
Two chiral centers → two dicts in the neighbours list.
"""

from pychemprojections.fischer.multiplechiral import (
    get_carbon_neighbours_info_multiple_chiral,
    get_right_and_left_substituents,
    remove_unnecessary_neighbours,
    get_condensed_form_info_of_substituents_multiple_chiral,
    get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain,
)

# Prepared SMILES for (2R,3S)-2,3-dibromobutane
PREPARED_SMILES = "[H]C([H])([H])[C@]([H])(Br)[C@@]([H])(Br)C([H])([H])[H]"


def _fresh_info():
    """Return a fresh copy of the carbon neighbours info (avoids mutation side-effects)."""
    return get_carbon_neighbours_info_multiple_chiral(PREPARED_SMILES)


class TestGetCarbonNeighboursInfoMultipleChiral:
    def test_returns_list(self):
        info = _fresh_info()
        assert isinstance(info, list)

    def test_two_chiral_centers_detected(self):
        info = _fresh_info()
        chiral = [c for c in info if c["carbon_type"] == "chiral"]
        assert len(chiral) == 2

    def test_each_entry_has_required_keys(self):
        info = _fresh_info()
        for entry in info:
            assert "carbon_type" in entry
            assert "carbon_atom_pos_in_string" in entry
            assert "substituents" in entry
            assert "config" in entry

    def test_first_chiral_center_has_four_substituents(self):
        info = _fresh_info()
        chiral_entries = [c for c in info if c["carbon_type"] == "chiral"]
        assert len(chiral_entries[0]["substituents"]) == 4

    def test_configurations_are_r_and_s(self):
        info = _fresh_info()
        configs = {c["config"] for c in info if c["carbon_type"] == "chiral"}
        assert configs == {"R", "S"}


class TestRemoveUnnecessaryNeighbours:
    def test_first_carbon_loses_last_substituent(self):
        info = _fresh_info()
        original_first_len = len(info[0]["substituents"])
        info = remove_unnecessary_neighbours(info)
        assert len(info[0]["substituents"]) == original_first_len - 1

    def test_last_carbon_loses_first_substituent(self):
        info = _fresh_info()
        original_last_len = len(info[-1]["substituents"])
        info = remove_unnecessary_neighbours(info)
        assert len(info[-1]["substituents"]) == original_last_len - 1

    def test_returns_same_number_of_entries(self):
        info = _fresh_info()
        original_len = len(info)
        info = remove_unnecessary_neighbours(info)
        assert len(info) == original_len


class TestGetRightAndLeftSubstituents:
    """
    After remove_unnecessary_neighbours the two chiral centers have 3 subs each.
    get_right_and_left_substituents must return two lists of equal length.
    """

    def _build_condensed_info(self):
        info = _fresh_info()
        info = remove_unnecessary_neighbours(info)
        info = get_condensed_form_info_of_substituents_multiple_chiral(info)
        return info

    def test_returns_two_lists(self):
        condensed = self._build_condensed_info()
        left, right = get_right_and_left_substituents(condensed)
        assert isinstance(left, list)
        assert isinstance(right, list)

    def test_left_and_right_have_same_length(self):
        condensed = self._build_condensed_info()
        left, right = get_right_and_left_substituents(condensed)
        assert len(left) == len(right)

    def test_substituents_are_strings(self):
        condensed = self._build_condensed_info()
        left, right = get_right_and_left_substituents(condensed)
        assert all(isinstance(s, str) for s in left)
        assert all(isinstance(s, str) for s in right)

    def test_length_equals_number_of_carbons(self):
        condensed = self._build_condensed_info()
        left, right = get_right_and_left_substituents(condensed)
        assert len(left) == len(condensed)


class TestGetSmilesOfChiralSubstituentGroupsInMultipleChiralChain:
    """Tests for the low-level bracket-slicing function."""

    def test_chiral_carbon_returns_four_groups(self):
        subs = get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
            PREPARED_SMILES,
            chiral_tag_of_c_atom=True,
            sq_bracket_begin_first_chiral_center=PREPARED_SMILES.index("[C@]"),
        )
        assert len(subs) == 4

    def test_returns_list_of_strings(self):
        subs = get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
            PREPARED_SMILES,
            chiral_tag_of_c_atom=True,
            sq_bracket_begin_first_chiral_center=PREPARED_SMILES.index("[C@]"),
        )
        assert all(isinstance(s, str) for s in subs)

    def test_at_H_variant_second_substituent_is_H(self):
        """When @H] is in the string, sub_group_b is hard-coded to 'H'."""
        smiles_with_H = "[C@H](Br)(Cl)F"
        subs = get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
            smiles_with_H,
            chiral_tag_of_c_atom=True,
            sq_bracket_begin_first_chiral_center=0,
        )
        assert subs[1] == "H"
