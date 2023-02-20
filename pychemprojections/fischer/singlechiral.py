from pychemprojections.utils.smiles_manipulation_utils import (
    get_index_of_corresponding_bracket,
)
from pychemprojections.fischer.stringmanipulation import (
    prepare_strings_for_fischer_projection_plot,
)
from typing import List
from pychemprojections.utils.logger_utils import get_module_logger
from pychemprojections.fischer.drawingclasses import DrawingInfo

logger = get_module_logger(__name__)


def get_smiles_of_chiral_substituent_groups_single_chiral(
    smiles_mol_prepared: str,
) -> List[str]:
    sq_bracket_begin = smiles_mol_prepared.find("[C@")
    sq_bracket_end = smiles_mol_prepared.find("@]") + 1

    sub_group_a = smiles_mol_prepared[:sq_bracket_begin]

    mol_b_normal_bracket_begin_idx = sq_bracket_end + 1
    mol_b_normal_bracket_end_idx = get_index_of_corresponding_bracket(
        smiles_mol_prepared, mol_b_normal_bracket_begin_idx
    )
    sub_group_b = smiles_mol_prepared[
        mol_b_normal_bracket_begin_idx + 1 : mol_b_normal_bracket_end_idx
    ]

    mol_c_normal_bracket_begin_idx = mol_b_normal_bracket_end_idx + 1

    mol_c_normal_bracket_end_idx = get_index_of_corresponding_bracket(
        smiles_mol_prepared, mol_c_normal_bracket_begin_idx
    )

    sub_group_c = smiles_mol_prepared[
        mol_c_normal_bracket_begin_idx + 1 : mol_c_normal_bracket_end_idx
    ]
    sub_group_d = smiles_mol_prepared[mol_c_normal_bracket_end_idx + 1 :]

    configuration = get_configuration_single_chiral_center(smiles_mol_prepared)
    logger.info(f"{configuration=}")
    if configuration == "S":
        substituent_groups = [sub_group_a, sub_group_b, sub_group_c, sub_group_d]
    else:
        substituent_groups = [sub_group_a, sub_group_c, sub_group_b, sub_group_d]

    return substituent_groups


def get_configuration_single_chiral_center(smiles_mol_prepared: str) -> str:
    if "@@" in smiles_mol_prepared:
        configuration = "S"
    else:
        configuration = "R"
    return configuration


def get_fisher_notation_of_all_substituents_single_chiral(
    substituents_condensed_form: List[str], drawing_info: DrawingInfo
) -> List[str]:
    return [
        prepare_strings_for_fischer_projection_plot(
            c, drawing_info.canvas_width_pixels, drawing_info.canvas_height_pixels
        )
        for c in substituents_condensed_form
    ]
