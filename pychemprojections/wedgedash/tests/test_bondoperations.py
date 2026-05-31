"""
Unit tests for pychemprojections.wedgedash.bondoperations

Tests cover:
  - get_reference_atom_details
  - get_bond_dir
  - set_wedge_bonds
  - collect_atom_ids_to_change_bond_type
  - add_missing_bond_dirs
"""

from rdkit import Chem
from rdkit.Chem import AllChem

from pychemprojections.wedgedash.bondoperations import (
    get_reference_atom_details,
    get_bond_dir,
    set_wedge_bonds,
    collect_atom_ids_to_change_bond_type,
    add_missing_bond_dirs,
    wedge,
    dash,
    in_plane_bond,
)
from pychemprojections.utils.rdkit_utils import (
    get_atom_ids_of_chiral_carbons,
    get_neighbours_of_all_chiral_C_atoms,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

CHIRAL_SMILES = "[C@@H](F)(Cl)Br"


def _wedged_mol(smiles: str):
    """Return a molecule with conformer and wedge bonds set."""
    mol = Chem.MolFromSmiles(smiles)
    return set_wedge_bonds(mol)


# ---------------------------------------------------------------------------
# set_wedge_bonds
# ---------------------------------------------------------------------------
class TestSetWedgeBonds:
    def test_returns_mol(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        assert mol is not None

    def test_molecule_has_conformer(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        assert mol.GetNumConformers() >= 1

    def test_hydrogens_added(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        h_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")
        assert h_count >= 1


# ---------------------------------------------------------------------------
# get_reference_atom_details
# ---------------------------------------------------------------------------
class TestGetReferenceAtomDetails:
    def test_returns_int(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        result = get_reference_atom_details(mol, chiral_ids[0])
        assert isinstance(result, int)

    def test_returns_valid_atom_idx(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        ref_idx = get_reference_atom_details(mol, chiral_ids[0])
        # Valid index: -1 means not found, otherwise within range
        if ref_idx != -1:
            assert 0 <= ref_idx < mol.GetNumAtoms()

    def test_returns_minus_one_when_no_wedge_bond(self):
        # Achiral molecule has no wedge bonds
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        result = get_reference_atom_details(mol, 0)
        assert result == -1


# ---------------------------------------------------------------------------
# get_bond_dir
# ---------------------------------------------------------------------------
class TestGetBondDir:
    def test_returns_bond_dir_type(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        result = get_bond_dir(mol, chiral_ids[0])
        assert result in (wedge, dash, in_plane_bond)

    def test_chiral_mol_has_non_in_plane_bond(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        result = get_bond_dir(mol, chiral_ids[0])
        assert result != in_plane_bond

    def test_achiral_carbon_returns_in_plane(self):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        result = get_bond_dir(mol, 0)
        assert result == in_plane_bond


# ---------------------------------------------------------------------------
# collect_atom_ids_to_change_bond_type
# ---------------------------------------------------------------------------
class TestCollectAtomIdsToChangeBondType:
    def test_returns_list(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        neigh_info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        result = collect_atom_ids_to_change_bond_type(mol, neigh_info)
        assert isinstance(result, list)

    def test_one_atom_per_chiral_center(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        neigh_info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        result = collect_atom_ids_to_change_bond_type(mol, neigh_info)
        assert len(result) == len(chiral_ids)

    def test_returned_ids_are_valid_atom_indices(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        neigh_info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        result = collect_atom_ids_to_change_bond_type(mol, neigh_info)
        for idx in result:
            assert 0 <= idx < mol.GetNumAtoms()


# ---------------------------------------------------------------------------
# add_missing_bond_dirs
# ---------------------------------------------------------------------------
class TestAddMissingBondDirs:
    def test_returns_mol(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        neigh_info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        atom_ids = collect_atom_ids_to_change_bond_type(mol, neigh_info)
        result = add_missing_bond_dirs(mol, atom_ids)
        assert result is not None

    def test_bond_direction_flipped_for_target_atoms(self):
        """The bond type at the target atom index should differ from in_plane_bond."""
        mol = _wedged_mol(CHIRAL_SMILES)
        chiral_ids = get_atom_ids_of_chiral_carbons(mol)
        neigh_info = get_neighbours_of_all_chiral_C_atoms(mol, chiral_ids)
        atom_ids = collect_atom_ids_to_change_bond_type(mol, neigh_info)
        result = add_missing_bond_dirs(mol, atom_ids)

        for bond in result.GetBonds():
            if bond.GetEndAtomIdx() in atom_ids:
                assert bond.GetBondDir() in (wedge, dash)

    def test_empty_ids_list_leaves_mol_unchanged(self):
        mol = _wedged_mol(CHIRAL_SMILES)
        # Record bond dirs before
        before = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondDir())
            for b in mol.GetBonds()
        ]
        add_missing_bond_dirs(mol, [])
        after = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondDir())
            for b in mol.GetBonds()
        ]
        assert before == after
