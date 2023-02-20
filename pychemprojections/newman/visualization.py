from pychemprojections.utils.rdkit_utils import (
    preprocess_molecule,
    get_iupac_name_from_smiles,
)
import os
import matplotlib.pyplot as plt
from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.newman.stringmanipulation import (
    prepare_strings_for_newman_projection_plot,
    get_condensed_groups,
    get_post_processed_smiles,
    create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def,
    get_substituents_front_and_rear,
    calibration_for_post_processing,
)
from pychemprojections.newman.drawingelements import draw_circle, add_atoms, add_lines

logger = get_module_logger(__name__)


def plot_projection(
    rear_atoms,
    front_atoms,
    rear_angles,
    front_angles,
    circle_radius,
    canvas_width_inches,
    canvas_height_inches,
    axes_x_limits,
    axes_y_limits,
    origin_x,
    origin_y,
    dpi,
    start_line_offset,
    end_line_offset,
    annotation_offset,
    iupac_name,
):
    output_img_dir = "../output_images"
    os.makedirs(output_img_dir, exist_ok=True)

    plt.ioff()
    fig, ax = plt.subplots(figsize=(canvas_width_inches, canvas_height_inches))
    fontsize = 20
    ax.set_xlim(*axes_x_limits)
    ax.set_ylim(*axes_y_limits)
    ax.axis("off")

    center = (origin_x, origin_y)
    ax = draw_circle(ax, center, circle_radius)
    logger.info("adding lines")
    pos = "rear"
    ax = add_lines(
        ax, rear_angles, pos, origin_x, origin_y, start_line_offset, end_line_offset
    )
    pos = "front"
    ax = add_lines(
        ax, front_angles, pos, origin_x, origin_y, start_line_offset, end_line_offset
    )
    logger.info("adding atoms")

    pos = "rear"
    ax = add_atoms(
        ax,
        rear_atoms,
        rear_angles,
        fontsize,
        origin_x,
        origin_y,
        annotation_offset,
        pos,
    )

    pos = "front"
    ax = add_atoms(
        ax,
        front_atoms,
        front_angles,
        fontsize,
        origin_x,
        origin_y,
        annotation_offset,
        pos,
    )

    if iupac_name:
        output_file_name = f"{iupac_name}_newman_projection.png"
    else:
        output_file_name = "output_newman.png"

    output_file_path = os.path.join(output_img_dir, output_file_name)

    plt.savefig(output_file_path, dpi=dpi)
    plt.close()


def plot_newman_projection(
    input_smiles,
    canvas_width_pixels,
    canvas_height_pixels,
    carbon_ids_bond_to_examine=None,
):
    origin_x = 0
    origin_y = 0
    end_line_offset = 1.4
    start_line_offset = 1
    n_carbons_for_truncation = 6
    annotation_offset = end_line_offset + 0.3
    radius = 1
    front_angles = [60, 180, 300]
    rear_angles = [0, 120, 240]
    dpi = 300
    canvas_width_inches = canvas_width_pixels / dpi
    canvas_height_inches = canvas_height_pixels / dpi

    axes_x_limits = (-end_line_offset * 2, end_line_offset * 2)
    axes_y_limits = (-end_line_offset * 2, end_line_offset * 2)

    smiles_mol_prepared, mol = preprocess_molecule(input_smiles)
    carbon_ids_mol_to_str_map = (
        create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
            mol, smiles_mol_prepared
        )
    )
    iupac_name = get_iupac_name_from_smiles(input_smiles)

    if carbon_ids_bond_to_examine == None:
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
    plot_projection(
        rear_atoms,
        front_atoms,
        rear_angles,
        front_angles,
        radius,
        canvas_width_inches,
        canvas_height_inches,
        axes_x_limits,
        axes_y_limits,
        origin_x,
        origin_y,
        dpi,
        start_line_offset,
        end_line_offset,
        annotation_offset,
        iupac_name,
    )


def get_front_and_back_groups_for_plotting(post_processed_smiles):
    newman_smiles = [
        prepare_strings_for_newman_projection_plot(psm) for psm in post_processed_smiles
    ]
    front_atoms = newman_smiles[:3]
    rear_atoms = newman_smiles[3:]

    return front_atoms, rear_atoms
