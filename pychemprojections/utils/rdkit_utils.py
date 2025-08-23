from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Mol
from rdkit.Chem import rdMolDescriptors
import pubchempy
import re
from typing import List, Tuple, Dict, Any
from pychemprojections.utils.logger_utils import get_module_logger

logger = get_module_logger(__name__)


def smiles_to_svg(
    molecule: Mol,
    svg_file_path: str,
    canvas_width_pixels: int,
    canvas_height_pixels: int,
):
    drawer = rdMolDraw2D.MolDraw2DSVG(canvas_width_pixels, canvas_height_pixels)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    with open(svg_file_path, "w") as f:
        f.write(svg)
    return svg


def convert_smiles_to_molecular_formula(input_smiles: str) -> str:
    mol = Chem.MolFromSmiles(input_smiles)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    return molecular_formula


def get_iupac_name_from_smiles(input_smiles: str) -> str:
    try:
        compounds = pubchempy.get_compounds(input_smiles, namespace="smiles")
        match = compounds[0]
        iupac_name = match.iupac_name
    except Exception:
        iupac_name = None
    return iupac_name


def smiles_to_condensed_form(input_smiles: str) -> str:
    regex = r"([Cc]|])+[0-9]"  # regex for removing cyclic numbering, it could c/C followed by a number or like [C@]/[C@@] followed by a number
    pos_of_cyclic_Cs = [
        substr.start() + 1 for substr in re.finditer(regex, input_smiles)
    ]
    input_smiles_list = list(input_smiles)

    for index in sorted(pos_of_cyclic_Cs, reverse=True):
        del input_smiles_list[index]

    smiles_removed_cyclic_subscipts = "".join(input_smiles_list)
    condensed_string = smiles_removed_cyclic_subscipts.replace("[H]", "H")
    condensed_string = condensed_string.replace("(H)", "H")
    condensed_string = condensed_string.replace("HHHH", "H4")
    condensed_string = condensed_string.replace("HHH", "H3")
    condensed_string = condensed_string.replace("HH", "H2")
    condensed_string = condensed_string.replace("-", "")
    condensed_string = condensed_string.replace("c", "C")
    condensed_string = condensed_string.replace("[C@]", "C")
    condensed_string = condensed_string.replace("[C@@]", "C")

    return condensed_string


def get_atom_ids_of_all_carbon_atoms(molecule: Mol) -> List[int]:
    atom_ids = []
    for a in molecule.GetAtoms():
        atom_id = a.GetIdx()
        atom_symbol = a.GetSymbol()
        if atom_symbol == "C":
            atom_ids.append(atom_id)
    return sorted(atom_ids)


def preprocess_molecule(input_smiles: str) -> Tuple[str, Mol]:
    mol = Chem.MolFromSmiles(input_smiles)
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    mol = Chem.AddHs(mol)
    smiles_mol_prepared = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles_mol_prepared)

    return smiles_mol_prepared, mol


def cleanup_Hs(input_smiles: str) -> str:
    if input_smiles[0:4] == "[H]C":
        input_smiles_H_cleaned = "C[H]" + input_smiles[4:]
    else:
        input_smiles_H_cleaned = input_smiles
    return input_smiles_H_cleaned


def get_neighbours_of_all_C_atoms(molecule: Mol) -> Dict[str, Any]:
    mol_info = {}
    for a in molecule.GetAtoms():
        atom_id = a.GetIdx()
        atom_symbol = a.GetSymbol()

        if atom_symbol == "C":
            neighbours = [
                (n.GetIdx(), n.GetSymbol(), n.GetAtomicNum()) for n in a.GetNeighbors()
            ]
            mol_info[atom_id] = {"neighbours": neighbours, "symbol": atom_symbol}

    return mol_info


def get_neighbours_of_all_chiral_C_atoms(
    molecule: Mol, atom_ids: List[int]
) -> Dict[int, Any]:
    mol_info = {}
    for a in molecule.GetAtoms():
        atom_id = a.GetIdx()
        if atom_id in atom_ids:
            neigh_atom_ids = [n.GetIdx() for n in a.GetNeighbors()]
            neigh_atom_symbols = [n.GetSymbol() for n in a.GetNeighbors()]
            mol_info[atom_id] = {
                "neigh_atom_ids": neigh_atom_ids,
                "neigh_atom_symbols": neigh_atom_symbols,
            }
    return mol_info


def get_atom_ids_of_chiral_carbons(mol: Mol) -> List[int]:
    chiral_centers = Chem.FindMolChiralCenters(mol)
    atom_ids = [info[0] for info in chiral_centers]
    return atom_ids
