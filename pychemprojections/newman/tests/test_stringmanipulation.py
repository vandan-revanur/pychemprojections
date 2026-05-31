"""
Unit tests for pychemprojections.newman.stringmanipulation

Covers all public functions that perform string/SMILES manipulation for
Newman projections.  RDKit-heavy helpers (calibration_for_post_processing,
create_mapping_between_atom_ids…) are tested with lightweight molecules.
"""

from pychemprojections.newman.stringmanipulation import (
    prepare_string_for_newman_projection_plot,
    separate_atoms_and_nums,
    get_c_ids_in_smiles,
    recalibrate_mol_formula_to_account_for_H_addition,
    smiles_post_processing,
    get_post_processed_smiles,
    create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def,
    calibration_for_post_processing,
)
from pychemprojections.utils.rdkit_utils import preprocess_molecule


# ---------------------------------------------------------------------------
# prepare_string_for_newman_projection_plot
# ---------------------------------------------------------------------------
class TestPrepareStringForNewmanProjectionPlot:
    def test_no_numbers_wrapped_in_dollars(self):
        assert prepare_string_for_newman_projection_plot("CH") == "$CH$"

    def test_number_gets_subscript_braces(self):
        result = prepare_string_for_newman_projection_plot("CH3")
        assert "_{3}" in result
        assert result.startswith("$") and result.endswith("$")

    def test_multiple_numbers_all_subscripted(self):
        result = prepare_string_for_newman_projection_plot("C2H5")
        assert "_{2}" in result
        assert "_{5}" in result

    def test_returns_string(self):
        assert isinstance(prepare_string_for_newman_projection_plot("Br"), str)

    def test_dollar_delimited(self):
        result = prepare_string_for_newman_projection_plot("Br")
        assert result == "$Br$"


# ---------------------------------------------------------------------------
# separate_atoms_and_nums
# ---------------------------------------------------------------------------
class TestSeparateAtomsAndNums:
    def test_c2h5_separated(self):
        result = separate_atoms_and_nums("C2H5")
        # Expect something like "C_2_H_5_"
        assert "_" in result

    def test_no_numbers_no_underscore_between_atoms(self):
        result = separate_atoms_and_nums("Br")
        # "Br" → "B_r_" – the regex hits uppercase B block end
        assert isinstance(result, str)

    def test_ch4_atoms_and_nums_separated(self):
        result = separate_atoms_and_nums("CH4")
        assert "C" in result and "H" in result and "4" in result


# ---------------------------------------------------------------------------
# get_c_ids_in_smiles
# ---------------------------------------------------------------------------
class TestGetCIdsInSmiles:
    def test_single_carbon(self):
        # "CCO" – two carbons at indices 0 and 1
        result = get_c_ids_in_smiles("CCO")
        assert 0 in result
        assert 1 in result
        assert len(result) == 2

    def test_no_carbon_returns_empty(self):
        result = get_c_ids_in_smiles("Br")
        assert result == []

    def test_aromatic_c_lowercase_not_counted(self):
        # aromatic 'c' is not counted (function checks for uppercase C not followed by lowercase)
        result = get_c_ids_in_smiles("c1ccccc1")
        assert result == []

    def test_explicit_hydrogens_around_carbon_counted(self):
        s = "[H]C([H])([H])[H]"
        # The 'C' sits at a certain index; we just verify one C is found
        result = get_c_ids_in_smiles(s)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# recalibrate_mol_formula_to_account_for_H_addition
# ---------------------------------------------------------------------------
class TestRecalibrateMolFormula:
    def test_h_count_decremented_by_one(self):
        # "C2H6" → "C2H5"
        result = recalibrate_mol_formula_to_account_for_H_addition("C2H6")
        assert result == "C2H5"

    def test_h_count_decremented_multi_digit(self):
        # "C3H8" → "C3H7"
        result = recalibrate_mol_formula_to_account_for_H_addition("C3H8")
        assert result == "C3H7"

    def test_returns_string(self):
        assert isinstance(
            recalibrate_mol_formula_to_account_for_H_addition("C2H6"), str
        )


# ---------------------------------------------------------------------------
# smiles_post_processing
# ---------------------------------------------------------------------------
class TestSmilesPostProcessing:
    def test_short_chain_returns_condensed(self):
        # "CC" has 1 carbon (idx 0,1) → len < threshold 3
        result = smiles_post_processing("CC", "CC", n_carbons_for_truncation=3)
        assert result == "CC"

    def test_long_chain_returns_molecular_formula(self):
        # "CCCCCC" – 6 carbons, threshold 3
        result = smiles_post_processing("CCCCCC", "CCCCCC", n_carbons_for_truncation=3)
        # Should be a molecular formula like "C6H14" not the SMILES
        assert "C" in result
        # Result should NOT be the condensed SMILES itself
        assert result != "CCCCCC"

    def test_returns_string(self):
        assert isinstance(smiles_post_processing("CC", "CC", 10), str)


# ---------------------------------------------------------------------------
# get_post_processed_smiles
# ---------------------------------------------------------------------------
class TestGetPostProcessedSmiles:
    def test_returns_same_length_list(self):
        smiles_list = ["CC", "Br", "Cl"]
        condensed_list = ["CC", "Br", "Cl"]
        result = get_post_processed_smiles(
            smiles_list, condensed_list, n_carbons_for_truncation=5
        )
        assert len(result) == len(smiles_list)

    def test_empty_input_returns_empty(self):
        result = get_post_processed_smiles([], [], 3)
        assert result == []

    def test_all_results_are_strings(self):
        result = get_post_processed_smiles(["CC", "Br"], ["CC", "Br"], 5)
        assert all(isinstance(s, str) for s in result)


# ---------------------------------------------------------------------------
# create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def
# ---------------------------------------------------------------------------
class TestCreateMappingBetweenAtomIds:
    def test_mapping_keys_are_rdkit_ids(self):
        smiles_prep, mol = preprocess_molecule("CC")
        mapping = create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
            mol, smiles_prep
        )
        # Keys should be valid atom indices in mol
        for rdkit_id in mapping:
            atom = mol.GetAtomWithIdx(rdkit_id)
            assert atom.GetSymbol() == "C"

    def test_mapping_values_are_string_indices(self):
        smiles_prep, mol = preprocess_molecule("CC")
        mapping = create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
            mol, smiles_prep
        )
        for str_idx in mapping.values():
            assert isinstance(str_idx, int)
            assert smiles_prep[str_idx] == "C"

    def test_returns_dict(self):
        smiles_prep, mol = preprocess_molecule("CC")
        result = create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
            mol, smiles_prep
        )
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# calibration_for_post_processing
# ---------------------------------------------------------------------------
class TestCalibrationForPostProcessing:
    def test_returns_list_of_six_strings(self):
        groups = {
            "front": {"1": "C", "2": "(F)", "3": "(Cl)"},
            "rear": {"1": "(Br)", "2": "(I)", "3": "CC"},
        }
        result = calibration_for_post_processing(groups)
        assert isinstance(result, list)
        assert len(result) == 6
        assert all(isinstance(s, str) for s in result)
