"""
Unit tests for pychemprojections.fischer.singlechiral

Molecules used as fixtures:
  R-config  → "[H][C@](F)(Cl)Br"   (@@-free  → code says R)
  S-config  → "[H][C@@](F)(Cl)Br"  (@@ present → code says S)

Note: the code resolves chirality via presence of "@@" in the string, not
through RDKit chiral perception, so the tests mirror that contract.
"""

from pychemprojections.fischer.singlechiral import (
    get_configuration_single_chiral_center,
    get_smiles_of_chiral_substituent_groups_single_chiral,
)

# Prepared SMILES produced by preprocess_molecule for our test molecules
SMILES_R = "[H][C@](F)(Cl)Br"  # no @@  → code returns "R"
SMILES_S = "[H][C@@](F)(Cl)Br"  # has @@ → code returns "S"


class TestGetConfigurationSingleChiralCenter:
    def test_double_at_sign_is_S(self):
        assert get_configuration_single_chiral_center(SMILES_S) == "S"

    def test_single_at_sign_is_R(self):
        assert get_configuration_single_chiral_center(SMILES_R) == "R"

    def test_molecule_without_at_sign_is_R(self):
        # Achiral SMILES – no '@' at all
        assert get_configuration_single_chiral_center("CCC") == "R"


class TestGetSmilesOfChiralSubstituentGroupsSingleChiral:
    """
    For "[H][C@@](F)(Cl)Br"  (S config):
      sub_a = "[H]"
      sub_b = "F"
      sub_c = "Cl"
      sub_d = "Br"
    S config → order kept as [a, b, c, d]
    """

    def test_s_config_returns_four_substituents(self):
        result = get_smiles_of_chiral_substituent_groups_single_chiral(SMILES_S)
        assert len(result) == 4

    def test_s_config_substituents_correct(self):
        result = get_smiles_of_chiral_substituent_groups_single_chiral(SMILES_S)
        # S: [a, b, c, d] = ["[H]", "F", "Cl", "Br"]
        assert result == ["[H]", "F", "Cl", "Br"]

    def test_r_config_swaps_b_and_c(self):
        """R config: order is [a, c, b, d] – b and c are swapped."""
        result = get_smiles_of_chiral_substituent_groups_single_chiral(SMILES_R)
        # R: [a, c, b, d] = ["[H]", "Cl", "F", "Br"]
        assert result == ["[H]", "Cl", "F", "Br"]

    def test_r_config_returns_four_substituents(self):
        result = get_smiles_of_chiral_substituent_groups_single_chiral(SMILES_R)
        assert len(result) == 4

    def test_substituents_are_strings(self):
        result = get_smiles_of_chiral_substituent_groups_single_chiral(SMILES_S)
        assert all(isinstance(s, str) for s in result)
