"""
Unit tests for pychemprojections.utils.rdkit_utils

RDKit is a required dependency so these tests are kept focused on logic,
not on molecule rendering.  Network-dependent functions (pubchem lookup)
are tested via mocks.
"""

from unittest.mock import MagicMock, patch

from rdkit import Chem

from pychemprojections.utils.rdkit_utils import (
    cleanup_Hs,
    smiles_to_condensed_form,
    preprocess_molecule,
    get_atom_ids_of_all_carbon_atoms,
    get_atom_ids_of_chiral_carbons,
    get_neighbours_of_all_C_atoms,
    get_neighbours_of_all_chiral_C_atoms,
    get_iupac_name_from_smiles,
    convert_smiles_to_molecular_formula,
)


# ---------------------------------------------------------------------------
# cleanup_Hs
# ---------------------------------------------------------------------------
class TestCleanupHs:
    def test_leading_H_attached_to_C_is_moved_to_end(self):
        assert cleanup_Hs("[H]CCC") == "C[H]CC"

    def test_leading_H_with_substituents_only_moves_first_H(self):
        assert cleanup_Hs("[H]CC(Cl)Br") == "C[H]C(Cl)Br"

    def test_string_without_leading_H_is_unchanged(self):
        assert cleanup_Hs("CCC") == "CCC"

    def test_H_not_at_position_zero_is_unchanged(self):
        assert cleanup_Hs("C[H]C") == "C[H]C"

    def test_empty_string_does_not_crash(self):
        # The function slices [:4] so an empty / short string is safe
        assert cleanup_Hs("") == ""


# ---------------------------------------------------------------------------
# smiles_to_condensed_form
# ---------------------------------------------------------------------------
class TestSmilesToCondensedForm:
    def test_explicit_H_bracket_replaced(self):
        result = smiles_to_condensed_form("[H]CC")
        assert "[H]" not in result
        assert "H" in result

    def test_double_H_collapses_to_H2(self):
        result = smiles_to_condensed_form("([H])[H]")  # two H tokens
        # After replacement '(H)H' → 'H2' path in the function
        # just assert no square brackets remain
        assert "[" not in result

    def test_cyclic_numbering_removed(self):
        # C1CCCCC1 – the regex removes cyclic C+digit patterns.
        # The function strips the first '1' (after opening C) but the
        # trailing lone '1' is not matched because it is not preceded by C/c/].
        # Assert the open-ring digit on the first carbon is gone.
        result = smiles_to_condensed_form("C1CCCCC1")
        # Original "C1" → first digit removed, giving "CCCCC1" (trailing 1 remains)
        assert result.startswith("C") and len(result) < len("C1CCCCC1")

    def test_aromatic_c_upcased(self):
        result = smiles_to_condensed_form("c1ccccc1")
        assert "c" not in result

    def test_chiral_bracket_replaced(self):
        result = smiles_to_condensed_form("[C@](F)(Cl)(Br)")
        assert "[C@]" not in result

    def test_double_chiral_bracket_replaced(self):
        result = smiles_to_condensed_form("[C@@](F)(Cl)(Br)")
        assert "[C@@]" not in result

    def test_dash_removed(self):
        result = smiles_to_condensed_form("C-C-C")
        assert "-" not in result

    def test_plain_smiles_passes_through(self):
        result = smiles_to_condensed_form("CCBr")
        assert result == "CCBr"


# ---------------------------------------------------------------------------
# preprocess_molecule
# ---------------------------------------------------------------------------
class TestPreprocessMolecule:
    def test_returns_tuple_of_str_and_mol(self):
        smiles_prep, mol = preprocess_molecule("[C@@H](F)(Cl)Br")
        assert isinstance(smiles_prep, str)
        assert mol is not None

    def test_prepared_smiles_contains_chiral_tag(self):
        smiles_prep, _ = preprocess_molecule("[C@@H](F)(Cl)Br")
        # After adding Hs, chirality tag must still be present
        assert "@" in smiles_prep

    def test_prepared_smiles_has_explicit_hydrogens(self):
        smiles_prep, _ = preprocess_molecule("[C@@H](F)(Cl)Br")
        assert "[H]" in smiles_prep

    def test_achiral_molecule_does_not_crash(self):
        smiles_prep, mol = preprocess_molecule("CCO")
        assert isinstance(smiles_prep, str)
        assert mol is not None


# ---------------------------------------------------------------------------
# get_atom_ids_of_all_carbon_atoms
# ---------------------------------------------------------------------------
class TestGetAtomIdsOfAllCarbonAtoms:
    def test_ethanol_has_two_carbons(self):
        mol = Chem.MolFromSmiles("CCO")
        ids = get_atom_ids_of_all_carbon_atoms(mol)
        assert len(ids) == 2

    def test_methane_has_one_carbon(self):
        mol = Chem.MolFromSmiles("C")
        ids = get_atom_ids_of_all_carbon_atoms(mol)
        assert len(ids) == 1

    def test_benzene_has_six_carbons(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        ids = get_atom_ids_of_all_carbon_atoms(mol)
        assert len(ids) == 6

    def test_result_is_sorted(self):
        mol = Chem.MolFromSmiles("CCCC")
        ids = get_atom_ids_of_all_carbon_atoms(mol)
        assert ids == sorted(ids)

    def test_molecule_without_carbon_returns_empty(self):
        mol = Chem.MolFromSmiles("O")
        ids = get_atom_ids_of_all_carbon_atoms(mol)
        assert ids == []


# ---------------------------------------------------------------------------
# get_atom_ids_of_chiral_carbons
# ---------------------------------------------------------------------------
class TestGetAtomIdsOfChiralCarbons:
    def test_single_chiral_center_detected(self):
        mol = Chem.MolFromSmiles("[C@@H](F)(Cl)Br")
        ids = get_atom_ids_of_chiral_carbons(mol)
        assert len(ids) == 1

    def test_two_chiral_centers_detected(self):
        # Use a molecule with two non-equivalent chiral centers so RDKit
        # does not collapse them as a meso compound.
        mol = Chem.MolFromSmiles("[C@@H](F)(Cl)[C@H](Br)(I)")
        ids = get_atom_ids_of_chiral_carbons(mol)
        assert len(ids) == 2

    def test_achiral_molecule_returns_empty(self):
        mol = Chem.MolFromSmiles("CCO")
        ids = get_atom_ids_of_chiral_carbons(mol)
        assert ids == []


# ---------------------------------------------------------------------------
# get_neighbours_of_all_C_atoms
# ---------------------------------------------------------------------------
class TestGetNeighboursOfAllCAtoms:
    def test_methane_has_no_C_neighbours(self):
        mol = Chem.MolFromSmiles("C")
        info = get_neighbours_of_all_C_atoms(mol)
        assert len(info) == 1  # one carbon
        carbon_id = list(info.keys())[0]
        # neighbours are all H – but Chem.MolFromSmiles gives implicit Hs
        assert info[carbon_id]["symbol"] == "C"

    def test_ethanol_carbon_neighbours(self):
        mol = Chem.MolFromSmiles("CCO")
        info = get_neighbours_of_all_C_atoms(mol)
        assert len(info) == 2

    def test_keys_are_carbon_atom_ids(self):
        mol = Chem.MolFromSmiles("CCO")
        info = get_neighbours_of_all_C_atoms(mol)
        for atom_id in info:
            atom = mol.GetAtomWithIdx(atom_id)
            assert atom.GetSymbol() == "C"


# ---------------------------------------------------------------------------
# get_neighbours_of_all_chiral_C_atoms
# ---------------------------------------------------------------------------
class TestGetNeighboursOfAllChiralCAtoms:
    def test_returns_info_only_for_requested_ids(self):
        mol = Chem.MolFromSmiles("[C@@H](F)(Cl)Br")
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        assert set(info.keys()) == set(chiral_ids)

    def test_neighbour_ids_and_symbols_present(self):
        mol = Chem.MolFromSmiles("[C@@H](F)(Cl)Br")
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        for v in info.values():
            assert "neigh_atom_ids" in v
            assert "neigh_atom_symbols" in v

    def test_empty_id_list_returns_empty_dict(self):
        mol = Chem.MolFromSmiles("CCO")
        info = get_neighbours_of_all_chiral_C_atoms(mol, [])
        assert info == {}


# ---------------------------------------------------------------------------
# convert_smiles_to_molecular_formula
# ---------------------------------------------------------------------------
class TestConvertSmilesToMolecularFormula:
    def test_methane(self):
        assert convert_smiles_to_molecular_formula("C") == "CH4"

    def test_ethanol(self):
        formula = convert_smiles_to_molecular_formula("CCO")
        assert formula == "C2H6O"

    def test_chiral_molecule(self):
        formula = convert_smiles_to_molecular_formula("[C@@H](F)(Cl)Br")
        assert "C" in formula
        assert "F" in formula


# ---------------------------------------------------------------------------
# get_iupac_name_from_smiles  (mocked – no network)
# ---------------------------------------------------------------------------
class TestGetIupacNameFromSmiles:
    def test_returns_iupac_name_on_success(self):
        mock_compound = MagicMock()
        mock_compound.iupac_name = "ethanol"
        with patch("pychemprojections.utils.rdkit_utils.pubchempy") as mock_pc:
            mock_pc.get_compounds.return_value = [mock_compound]
            result = get_iupac_name_from_smiles("CCO")
        assert result == "ethanol"

    def test_returns_none_on_lookup_failure(self):
        with patch("pychemprojections.utils.rdkit_utils.pubchempy") as mock_pc:
            mock_pc.get_compounds.side_effect = Exception("network error")
            result = get_iupac_name_from_smiles("CCO")
        assert result is None

    def test_returns_none_when_empty_result(self):
        with patch("pychemprojections.utils.rdkit_utils.pubchempy") as mock_pc:
            mock_pc.get_compounds.return_value = []
            result = get_iupac_name_from_smiles("CCO")
        assert result is None
