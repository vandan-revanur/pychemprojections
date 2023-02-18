import matplotlib.pyplot as plt
import os
from rdkit import Chem
from pychemprojections.utils.rdkit_utils import (
    get_iupac_name_from_smiles,
    preprocess_molecule,
)
from pychemprojections.fischer.singlechiral import (
    get_smiles_of_chiral_substituent_groups_single_chiral,
    get_fisher_notation_of_all_substituents_single_chiral,
)
from pychemprojections.fischer.stringmanipulation import (
    get_condensed_smiles_of_all_substituents,
)
from pychemprojections.fischer.multiplechiral import (
    get_carbon_neighbours_info_multiple_chiral,
    get_right_and_left_substituents,
    get_condensed_form_info_of_substituents_multiple_chiral,
    remove_unnecessary_neighbours,
    get_fisher_notation_info_for_substituents_multiple_chiral,
)
from pychemprojections.utils.rdkit_utils import cleanup_Hs
from typing import List
from pychemprojections.utils.logger_utils import get_module_logger

logger = get_module_logger(__name__)


def plot_fisher_projection_multiple_chiral_centers(
    up: str,
    down: str,
    left: List[str],
    right: List[str],
    iupac_name: str,
    canvas_width: int,
    canvas_height: int,
    output_img_path: str = "output.png",
):
    n_right = len(right)
    right = reversed(right)
    left = reversed(left)

    DPI = 300
    line_color = "black"
    line_1_vertical_x_coords = (canvas_width / 2, canvas_width / 2)
    line_1_vertical_y_coords = (0, canvas_height)

    fig, ax1 = plt.subplots(figsize=(canvas_width / DPI, canvas_height / DPI))

    ytick_locations = []
    line_2_horizontal_x_coords = (0, canvas_width)

    for i in range(1, n_right + 1):
        line_2_horizontal_y_coords = (
            canvas_height * (i) / (n_right + 1),
            canvas_height * (i) / (n_right + 1),
        )
        ytick_locations.append(canvas_height * (i) / (n_right + 1))
        ax1.plot(
            line_1_vertical_x_coords,
            line_1_vertical_y_coords,
            line_2_horizontal_x_coords,
            line_2_horizontal_y_coords,
            c=line_color,
        )

    ax1.set_ylim(0, canvas_height)
    ax1.set_xlim(0, canvas_width)

    ax1.set_xticks((canvas_width / 2,))
    ax1.set_yticks(ytick_locations)

    ax1.tick_params(color="white")
    labels_left = ax1.set_yticklabels(left)
    labels_down = ax1.set_xticklabels((down,), ha="left")

    # y-axis on right
    ax2 = ax1.twinx()
    ax2.tick_params(color="white")
    ax2.set_ylim(0, canvas_height)
    ax2.set_xlim(0, canvas_width)
    ax2.set_yticks(ytick_locations)
    labels_right = ax2.set_yticklabels(right)

    # x-axis on top
    ax3 = ax1.twiny()

    ax3.tick_params(color="white")
    ax3.set_ylim(0, canvas_height)
    ax3.set_xlim(0, canvas_width)
    ax3.set_xticks((canvas_width / 2,))
    labels_up = ax3.set_xticklabels((up,), ha="left", wrap=True)
    if iupac_name:
        plt.title(iupac_name + "\n")

    for key in ax1.spines.keys():
        ax1.spines[key].set_visible(False)

    for key in ax2.spines.keys():
        ax2.spines[key].set_visible(False)

    for key in ax3.spines.keys():
        ax3.spines[key].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_img_path)


def plot_fisher_projection_single_chiral_center(
    up: str,
    down: str,
    left: str,
    right: str,
    iupac_name: str,
    canvas_width: int,
    canvas_height: int,
    output_img_path: str = "output_single_chiral.png",
):
    DPI = 300
    line_color = "black"
    line_1_vertical_x_coords = (canvas_width / 2, canvas_width / 2)
    line_1_vertical_y_coords = (0, canvas_height)

    line_2_horizontal_x_coords = (0, canvas_width)
    line_2_horizontal_y_coords = (canvas_height / 2, canvas_height / 2)
    fig, ax1 = plt.subplots(figsize=(canvas_width / DPI, canvas_height / DPI))

    ax1.plot(
        line_1_vertical_x_coords,
        line_1_vertical_y_coords,
        line_2_horizontal_x_coords,
        line_2_horizontal_y_coords,
        c=line_color,
    )
    ax1.set_ylim(0, canvas_height)
    ax1.set_xlim(0, canvas_width)

    ax1.set_yticks((canvas_height / 2,))
    ax1.set_xticks((canvas_width / 2,))
    ax1.tick_params(color="white")
    labels_left = ax1.set_yticklabels((left,))
    labels_down = ax1.set_xticklabels((down,), ha="left")

    # x-axis on right
    ax2 = ax1.twinx()

    ax2.tick_params(color="white")
    ax2.set_ylim(0, canvas_height)
    ax2.set_xlim(0, canvas_width)
    ax2.set_yticks((canvas_height / 2,))
    labels_right = ax2.set_yticklabels((right,))

    # x-axis on top
    ax3 = ax1.twiny()

    ax3.tick_params(color="white")
    ax3.set_ylim(0, canvas_height)
    ax3.set_xlim(0, canvas_width)
    ax3.set_xticks((canvas_width / 2,))
    labels_up = ax3.set_xticklabels((up,), ha="left", wrap=True)

    if iupac_name:
        plt.title(iupac_name + "\n")

    for key in ax1.spines.keys():
        ax1.spines[key].set_visible(False)

    for key in ax2.spines.keys():
        ax2.spines[key].set_visible(False)

    for key in ax3.spines.keys():
        ax3.spines[key].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_img_path)


def plot_fischer_projection(
    input_smiles: str, canvas_width: int = 1000, canvas_height: int = 1000
):
    output_img_dir = "../output_images"
    os.makedirs(output_img_dir, exist_ok=True)

    iupac_name = get_iupac_name_from_smiles(input_smiles)
    smiles_mol_prepared, mol = preprocess_molecule(input_smiles)

    logger.info("Preparing SMILES for drawing")
    logger.debug(f"{smiles_mol_prepared=}")
    chiral_centers = Chem.FindMolChiralCenters(mol)
    n_chiral_centers = len(chiral_centers)
    logger.info("num of chiral centers %s", n_chiral_centers)

    if n_chiral_centers == 1:
        logger.info("Extracting the SMILES of the substituents")
        substituent_groups = get_smiles_of_chiral_substituent_groups_single_chiral(
            smiles_mol_prepared
        )
        logger.debug(f"{substituent_groups=}")

        logger.info("Converting the SMILES substituents into condensed form")
        substituents_condensed_form = get_condensed_smiles_of_all_substituents(
            substituent_groups
        )
        logger.debug(f"{substituents_condensed_form=}")

        logger.info(
            "Converting the condensed notation into a format suitable for plotting on fischer diagram"
        )
        fischer_notation = get_fisher_notation_of_all_substituents_single_chiral(
            substituents_condensed_form, canvas_width, canvas_height
        )
        logger.debug(f"{fischer_notation=}")

        [up, right, left, down] = fischer_notation
        if iupac_name:
            output_file_name = f"{iupac_name}_fischer_projection.png"
        else:
            output_file_name = "output_single_chiral_center.png"

        output_file_path = os.path.join(output_img_dir, output_file_name)
        logger.info("Plotting the fischer projection diagram")
        plot_fisher_projection_single_chiral_center(
            up=up,
            down=down,
            left=left,
            right=right,
            iupac_name=iupac_name,
            canvas_width=canvas_width,
            canvas_height=canvas_height,
            output_img_path=output_file_path,
        )
    else:
        smiles_mol_prepared = cleanup_Hs(smiles_mol_prepared)
        logger.info("Extracting the SMILES of the substituents")
        carbons_neighbours_info = get_carbon_neighbours_info_multiple_chiral(
            smiles_mol_prepared
        )
        logger.debug(f"{carbons_neighbours_info=}")

        sorted_carbons_neighbours_info = sorted(
            carbons_neighbours_info, key=lambda d: d["carbon_atom_pos_in_string"]
        )

        sorted_carbons_neighbours_info = remove_unnecessary_neighbours(
            sorted_carbons_neighbours_info
        )

        logger.info("Converting the SMILES substituents into condensed form")
        substituents_condensed_form = (
            get_condensed_form_info_of_substituents_multiple_chiral(
                sorted_carbons_neighbours_info
            )
        )
        logger.debug(f"{substituents_condensed_form=}")

        logger.info(
            "Converting the condensed notation into a format suitable for plotting on fischer diagram"
        )
        substituents_fischer_notation = (
            get_fisher_notation_info_for_substituents_multiple_chiral(
                substituents_condensed_form, canvas_width, canvas_height
            )
        )
        logger.debug(f"{substituents_fischer_notation=}")

        up = substituents_fischer_notation[0]["substituents"][0]
        down = substituents_fischer_notation[-1]["substituents"][1]
        left, right = get_right_and_left_substituents(substituents_condensed_form)
        assert len(left) == len(
            right
        ), f"The number of substituents on the left side of the chain ({len(left)}) do not match the number of substituents on the right side ({len(right)})"

        if iupac_name:
            output_file_name = f"{iupac_name}_fischer_projection.png"
        else:
            output_file_name = "output_multiple_chiral_centers.png"

        output_file_path = os.path.join(output_img_dir, output_file_name)

        logger.info("Plotting the fischer projection diagram")
        plot_fisher_projection_multiple_chiral_centers(
            up=up,
            down=down,
            left=left,
            right=right,
            iupac_name=iupac_name,
            canvas_width=canvas_width,
            canvas_height=canvas_height,
            output_img_path=output_file_path,
        )
