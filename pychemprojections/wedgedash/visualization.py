from pychemprojections.utils.rdkit_utils import (
    preprocess_molecule,
    get_iupac_name_from_smiles,
)
from rdkit.Chem import Draw
import os
from rdkit.Chem import AllChem
from pychemprojections.utils.rdkit_utils import smiles_to_svg
from pychemprojections.utils.rdkit_utils import (
    get_atom_ids_of_chiral_carbons,
    get_neighbours_of_all_chiral_C_atoms,
)
from pychemprojections.wedgedash.bondoperations import (
    collect_atom_ids_to_change_bond_type,
    add_missing_bond_dirs,
    set_wedge_bonds,
)


def plot_wedgedash_projection(
    input_smiles: str,
    canvas_width_pixels: int = 500,
    canvas_height_pixels: int = 500,
):
    """
    Plot/Draw the WedgeDash projection

    Parameters
    ----------
    input_smiles : Mol
        Molecule to examine

    canvas_width_pixels: int
        Width of the canvas in pixels

    canvas_height_pixels: int
        Height of the canvas in pixels

    """
    iupac_name = get_iupac_name_from_smiles(input_smiles)
    _, mol = preprocess_molecule(input_smiles)
    mol = set_wedge_bonds(mol)
    draw_opts = Draw.MolDrawOptions()
    draw_opts.addAtomIndices = True
    atom_ids = get_atom_ids_of_chiral_carbons(mol)
    substituent_neighbours_chiral_carbons = get_neighbours_of_all_chiral_C_atoms(
        mol, atom_ids
    )
    atoms_ids_to_change_bond_type = collect_atom_ids_to_change_bond_type(
        mol, substituent_neighbours_chiral_carbons
    )
    mol = add_missing_bond_dirs(mol, atoms_ids_to_change_bond_type)
    AllChem.EmbedMolecule(mol)

    output_img_dir = "output_images"
    os.makedirs(output_img_dir, exist_ok=True)
    if iupac_name:
        output_file_name = f"{iupac_name}_wedgedash_projection.svg"
    else:
        output_file_name = "output_wedgedash.svg"

    output_file_path = os.path.join(output_img_dir, output_file_name)

    output_img = smiles_to_svg(
        mol, output_file_path, canvas_width_pixels, canvas_height_pixels
    )
    return output_img
