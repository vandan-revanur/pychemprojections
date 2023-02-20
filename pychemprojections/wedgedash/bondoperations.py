from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Mol
from pychemprojections.utils.logger_utils import get_module_logger
from typing import Dict, Any, List

from rdkit.Chem import PropertyMol

logger = get_module_logger(__name__)

wedge = Chem.rdchem.BondDir.BEGINWEDGE
dash = Chem.rdchem.BondDir.BEGINDASH
in_plane_bond = Chem.rdchem.BondDir.NONE


def get_reference_atom_details(
    mol: PropertyMol, corresponding_chiral_carbon_idx: int
) -> int:
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtomIdx() == corresponding_chiral_carbon_idx
            and bond.GetBondDir() != in_plane_bond
        ):
            return bond.GetEndAtomIdx()
        else:
            pass
    return -1


def get_bond_dir(
    mol: PropertyMol, corresponding_chiral_carbon_idx: int
) -> Chem.rdchem.BondDir:
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtomIdx() == corresponding_chiral_carbon_idx
            and bond.GetBondDir() != in_plane_bond
        ):
            return bond.GetBondDir()
        else:
            pass
    return in_plane_bond


def set_wedge_bonds(mol: Mol) -> Mol:
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformers()[0])
    return mol


def collect_atom_ids_to_change_bond_type(
    mol: PropertyMol, substituent_neighbours_chiral_carbons: Dict[int, Any]
) -> List[int]:
    atoms_ids_to_change_bond_type = []
    for chiral_carbon_atom_idx, n_info in substituent_neighbours_chiral_carbons.items():
        reference_atom_idx = get_reference_atom_details(mol, chiral_carbon_atom_idx)
        neighbour_idx_in_group = n_info["neigh_atom_ids"].index(reference_atom_idx)

        if neighbour_idx_in_group == 1:
            atom_idx_in_neighbours = 2
        elif neighbour_idx_in_group == 2:
            atom_idx_in_neighbours = 1
        elif neighbour_idx_in_group == 3:
            atom_idx_in_neighbours = 1
        else:
            atom_idx_in_neighbours = 1

        atom_idx = n_info["neigh_atom_ids"][atom_idx_in_neighbours]
        atoms_ids_to_change_bond_type.append(atom_idx)

    return atoms_ids_to_change_bond_type


def add_missing_bond_dirs(
    mol: PropertyMol, atoms_ids_to_change_bond_type: List[int]
) -> PropertyMol:
    for bond in mol.GetBonds():
        end_atom_idx = bond.GetEndAtomIdx()

        if end_atom_idx in atoms_ids_to_change_bond_type:
            corresponding_chiral_carbon_idx = bond.GetBeginAtomIdx()
            reference_bond_type = get_bond_dir(mol, corresponding_chiral_carbon_idx)

            if reference_bond_type == dash:
                Chem.rdchem.Bond.SetBondDir(bond, wedge)
            elif reference_bond_type == wedge:
                Chem.rdchem.Bond.SetBondDir(bond, dash)
            else:
                logger.warning(
                    "reference_bond_type not recognized %s", reference_bond_type
                )

    return mol
