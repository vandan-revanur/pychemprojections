from pychemprojections.utils.rdkit_utils import (
    preprocess_molecule,
    get_iupac_name_from_smiles,
)
import os
import matplotlib.pyplot as plt
from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.newman.stringmanipulation import (
    prepare_string_for_newman_projection_plot,
    get_condensed_groups,
    get_post_processed_smiles,
    create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def,
    get_substituents_front_and_rear,
    calibration_for_post_processing,
)
from pychemprojections.newman.drawingelements import draw_circle, add_atoms, add_lines
from typing import Tuple, List
from pychemprojections.newman.drawingclasses import NewmanDrawingInfo, AtomInfo

logger = get_module_logger(__name__)


def plot_newman_projection(newman_drawing_info: NewmanDrawingInfo):
    """
    Plot the Newman projection

    Parameters
    ----------
    newman_drawing_info : NewmanDrawingInfo
        Infomation neccessary for plotting/drawing the newman projection such as front and rear atom labels, canvas
        width and height, iupac name of the substituents

    """
    output_img_dir = "output_images"
    origin_x = 0
    origin_y = 0
    end_line_offset = 1.4
    start_line_offset = 1

    annotation_offset = end_line_offset + 0.3
    front_angles = [60, 180, 300]
    rear_angles = [0, 120, 240]
    dpi = 300
    circle_radius = 1

    canvas_width_pixels = newman_drawing_info.canvas_width_pixels
    canvas_height_pixels = newman_drawing_info.canvas_height_pixels
    rear_atoms = newman_drawing_info.rear_atoms
    front_atoms = newman_drawing_info.front_atoms
    iupac_name = newman_drawing_info.iupac_name

    canvas_width_inches = canvas_width_pixels / dpi
    canvas_height_inches = canvas_height_pixels / dpi

    axes_x_limits = (-end_line_offset * 2, end_line_offset * 2)
    axes_y_limits = (-end_line_offset * 2, end_line_offset * 2)

    os.makedirs(output_img_dir, exist_ok=True)

    plt.ioff()
    fig, ax = plt.subplots(figsize=(canvas_width_inches, canvas_height_inches))
    fontsize = 20
    ax.set_xlim(*axes_x_limits)
    ax.set_ylim(*axes_y_limits)
    ax.axis("off")

    center = (origin_x, origin_y)
    ax = draw_circle(ax, center, circle_radius)

    rear_atom_info = AtomInfo(
        ax,
        rear_atoms,
        rear_angles,
        fontsize,
        origin_x,
        origin_y,
        annotation_offset,
        "rear",
    )

    front_atom_info = AtomInfo(
        ax,
        front_atoms,
        front_angles,
        fontsize,
        origin_x,
        origin_y,
        annotation_offset,
        "front",
    )
    logger.info("adding lines")
    ax = add_lines(rear_atom_info, start_line_offset, end_line_offset)
    ax = add_lines(front_atom_info, start_line_offset, end_line_offset)

    logger.info("adding atoms")
    ax = add_atoms(front_atom_info)
    ax = add_atoms(rear_atom_info)

    if iupac_name:
        output_file_name = f"{iupac_name}_newman_projection.png"
    else:
        output_file_name = "output_newman.png"

    output_file_path = os.path.join(output_img_dir, output_file_name)

    plt.savefig(output_file_path, dpi=dpi)
    plt.show()


def get_newman_drawing_info(
    input_smiles: str,
    canvas_width_pixels: int = 2000,
    canvas_height_pixels: int = 2000,
    carbon_ids_bond_to_examine: Tuple[int, int] = None,
) -> NewmanDrawingInfo:
    """
    Get all the information required to make a Newman Projection plot

    Parameters
    ----------
    input_smiles : str
        Input SMILES of the compound to examine
    canvas_width_pixels: int
        Width of the canvas to draw in terms of pixels
    canvas_height_pixels: int
        Height of the canvas to draw in terms of pixels
    carbon_ids_bond_to_examine:Tuple[int, int]
        Ids of carbon atoms at which we are observing the molecule in order to make the Newman projection

    Returns
    -------
    NewmanDrawingInfo
    Infomation neccessary for plotting/drawing the newman projection such as front and rear atom labels, canvas
        width and height, iupac name of the substituents

    """
    n_carbons_for_truncation = 6

    smiles_mol_prepared, mol = preprocess_molecule(input_smiles)
    carbon_ids_mol_to_str_map = (
        create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
            mol, smiles_mol_prepared
        )
    )
    iupac_name = get_iupac_name_from_smiles(input_smiles)

    if carbon_ids_bond_to_examine is None:
        start_carbon_mol_id_for_examination = list(carbon_ids_mol_to_str_map.keys())[0]
        end_carbon_mol_id_for_examination = list(carbon_ids_mol_to_str_map.keys())[1]

    else:
        start_carbon_mol_id_for_examination = carbon_ids_bond_to_examine[0]
        end_carbon_mol_id_for_examination = carbon_ids_bond_to_examine[1]

    start_carbon_atom_idx_in_str = carbon_ids_mol_to_str_map[
        start_carbon_mol_id_for_examination
    ]
    end_carbon_atom_idx_in_str = carbon_ids_mol_to_str_map[
        end_carbon_mol_id_for_examination
    ]
    groups = get_substituents_front_and_rear(
        smiles_mol_prepared, start_carbon_atom_idx_in_str, end_carbon_atom_idx_in_str
    )
    groups_condensed = get_condensed_groups(groups)
    smiles_groups = calibration_for_post_processing(groups)
    post_processed_smiles = get_post_processed_smiles(
        smiles_groups, groups_condensed, n_carbons_for_truncation
    )
    front_atoms, rear_atoms = get_front_and_back_groups_for_plotting(
        post_processed_smiles
    )

    newman_drawing_info = NewmanDrawingInfo(
        rear_atoms,
        front_atoms,
        canvas_width_pixels,
        canvas_height_pixels,
        iupac_name,
    )
    return newman_drawing_info


def get_front_and_back_groups_for_plotting(
    post_processed_smiles: List[str],
) -> Tuple[List[str], List[str]]:
    """
    Get the front and the back groups of substituents for plotting

    Parameters
    ----------
    post_processed_smiles : type
        SMILES post processed after replacing with molecular formula where appropriate

    Returns
    -------
    Tuple[List[str], List[str]]
    Tuple with a list of the front and list the back atoms

    """
    newman_smiles = [
        prepare_string_for_newman_projection_plot(psm) for psm in post_processed_smiles
    ]
    front_atoms = newman_smiles[:3]
    rear_atoms = newman_smiles[3:]

    return front_atoms, rear_atoms
