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

from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.fischer.drawingclasses import (
    MultipleChiralFischerNotation,
    DrawingInfo,
    SingleChiralFischerNotation,
)


logger = get_module_logger(__name__)


def plot_fisher_projection_multiple_chiral_centers(
    multiple_chiral_fischer_notation: MultipleChiralFischerNotation,
    drawing_info: DrawingInfo,
):
    """
    Plot/Draw the Fischer projection of the compound with multiple chiral carbons

    Parameters
    ----------
    multiple_chiral_fischer_notation : MultipleChiralFischerNotation
        MultipleChiralFischerNotation class contains the substituents strings on the up, down, left and right side of the
        Fischer projection backbone/main chain

    drawing_info: DrawingInfo
        DrawingInfo class contains information useful for drawing the Fischer projection such as
        input_smiles,canvas_width_pixels, canvas_height_pixels,smiles_preprocessed, and iupac_name

    """
    right = multiple_chiral_fischer_notation.right
    left = multiple_chiral_fischer_notation.left
    down = multiple_chiral_fischer_notation.down
    up = multiple_chiral_fischer_notation.up

    canvas_width = drawing_info.canvas_width_pixels
    canvas_height = drawing_info.canvas_height_pixels
    iupac_name = drawing_info.iupac_name

    output_img_dir = "../output_images"
    os.makedirs(output_img_dir, exist_ok=True)

    if iupac_name:
        output_file_name = f"{iupac_name}_fischer_projection.png"
    else:
        output_file_name = "output_single_chiral_center.png"

    output_img_path = os.path.join(output_img_dir, output_file_name)

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
    ax1.set_yticklabels(left)
    ax1.set_xticklabels((down,), ha="left")

    # y-axis on right
    ax2 = ax1.twinx()
    ax2.tick_params(color="white")
    ax2.set_ylim(0, canvas_height)
    ax2.set_xlim(0, canvas_width)
    ax2.set_yticks(ytick_locations)
    ax2.set_yticklabels(right)

    # x-axis on top
    ax3 = ax1.twiny()

    ax3.tick_params(color="white")
    ax3.set_ylim(0, canvas_height)
    ax3.set_xlim(0, canvas_width)
    ax3.set_xticks((canvas_width / 2,))
    ax3.set_xticklabels((up,), ha="left", wrap=True)
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
    single_chiral_fischer_notation: SingleChiralFischerNotation,
    drawing_info: DrawingInfo,
):
    """
    Plot/Draw the Fischer projection of the compound with single chiral carbon

    Parameters
    ----------
    single_chiral_fischer_notation : SingleChiralFischerNotation
        SingleChiralFischerNotation class contains the substituents strings on the up, down, left and right side of the
        chiral carbon

    drawing_info: DrawingInfo
        DrawingInfo class contains information useful for drawing the Fischer projection such as
        input_smiles,canvas_width_pixels, canvas_height_pixels,smiles_preprocessed, and iupac_name

    """
    DPI = 300
    line_color = "black"
    canvas_width = drawing_info.canvas_width_pixels
    canvas_height = drawing_info.canvas_height_pixels
    iupac_name = drawing_info.iupac_name

    up = single_chiral_fischer_notation.up
    down = single_chiral_fischer_notation.down
    left = single_chiral_fischer_notation.left
    right = single_chiral_fischer_notation.right

    output_img_dir = "pychemprojections/output_images"
    os.makedirs(output_img_dir, exist_ok=True)

    if iupac_name:
        output_file_name = f"{iupac_name}_fischer_projection.png"
    else:
        output_file_name = "output_single_chiral_center.png"

    output_img_path = os.path.join(output_img_dir, output_file_name)

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
    ax1.set_yticklabels((left,))
    ax1.set_xticklabels((down,), ha="left")

    # x-axis on right
    ax2 = ax1.twinx()

    ax2.tick_params(color="white")
    ax2.set_ylim(0, canvas_height)
    ax2.set_xlim(0, canvas_width)
    ax2.set_yticks((canvas_height / 2,))
    ax2.set_yticklabels((right,))

    # x-axis on top
    ax3 = ax1.twiny()

    ax3.tick_params(color="white")
    ax3.set_ylim(0, canvas_height)
    ax3.set_xlim(0, canvas_width)
    ax3.set_xticks((canvas_width / 2,))
    ax3.set_xticklabels((up,), ha="left", wrap=True)

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
    input_smiles: str,
    canvas_width: int = 1000,
    canvas_height: int = 1000,
):
    """
    XXXXXXX

    Parameters
    ----------
    input_smiles : str
        SMILES string of the compound whose Fischer Projection should be drawn

    canvas_width: int
        Width of the canvas to draw the Fischer Projection in pixels
    canvas_height: int
        Height of the canvas to draw the Fischer Projection in pixels.

    """
    iupac_name = get_iupac_name_from_smiles(input_smiles)
    smiles_mol_prepared, mol = preprocess_molecule(input_smiles)

    logger.info("Preparing SMILES for drawing")
    logger.debug(f"{smiles_mol_prepared=}")
    chiral_centers = Chem.FindMolChiralCenters(mol)
    n_chiral_centers = len(chiral_centers)
    logger.info("num of chiral centers %s", n_chiral_centers)

    drawing_info = DrawingInfo(
        input_smiles, canvas_width, canvas_height, smiles_mol_prepared, iupac_name
    )

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
            substituents_condensed_form, drawing_info
        )
        logger.debug(f"{fischer_notation=}")

        [up, right, left, down] = fischer_notation
        single_chiral_fischer_notation = SingleChiralFischerNotation(
            up=up, down=down, left=left, right=right
        )

        logger.info("Plotting the fischer projection diagram")
        plot_fisher_projection_single_chiral_center(
            single_chiral_fischer_notation, drawing_info
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
                substituents_condensed_form, drawing_info
            )
        )
        logger.debug(f"{substituents_fischer_notation=}")

        up = substituents_fischer_notation[0]["substituents"][0]
        down = substituents_fischer_notation[-1]["substituents"][1]
        left, right = get_right_and_left_substituents(substituents_condensed_form)
        assert len(left) == len(right), (
            f"The number of substituents on the left side of the chain ({len(left)}) do not match the number of substituents on the right side ({len(right)})"
        )

        multiple_chiral_fischer_notation = MultipleChiralFischerNotation(
            up=up, down=down, left=left, right=right
        )

        logger.info("Plotting the fischer projection diagram")
        plot_fisher_projection_multiple_chiral_centers(
            multiple_chiral_fischer_notation, drawing_info
        )
