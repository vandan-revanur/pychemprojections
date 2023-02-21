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
    """
    Get the SMILES of all the substituent groups attached to the chiral carbon atom

    Parameters
    ----------
    smiles_mol_prepared : str
        The SMILES of the preprocessed molecule.
        For details related to preprocessing see the function: pychemprojections.utils.rdkit_utils.preprocess_molecule

    Returns
    -------
    List[str]
    List of the four substituent groups attached to the chiral carbon atom.

    """
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
    if configuration == "S":
        substituent_groups = [sub_group_a, sub_group_b, sub_group_c, sub_group_d]
    else:
        substituent_groups = [sub_group_a, sub_group_c, sub_group_b, sub_group_d]

    return substituent_groups


def get_configuration_single_chiral_center(smiles_mol_prepared: str) -> str:
    """
    Get the chirality configuration at the chiral carbon

    Parameters
    ----------
    smiles_mol_prepared : str
        The SMILES of the preprocessed molecule.
        For details related to preprocessing see the function: pychemprojections.utils.rdkit_utils.preprocess_molecule

    Returns
    -------
    str
    Chirality configuration at the chiral configuration. Either R or S.

    """
    if "@@" in smiles_mol_prepared:
        configuration = "S"
    else:
        configuration = "R"
    return configuration


def get_fisher_notation_of_all_substituents_single_chiral(
    substituents_condensed_form: List[str], drawing_info: DrawingInfo
) -> List[str]:
    """
    Convert the condensed form of representation to the chemical notation with subscript of atom numbers for the
    Fischer projection plot

    Parameters
    ----------
    substituents_condensed_form : List[Dict[str, Any]]
        List of the information dictionaries that consist of the condensed form of the substituents
    drawing_info: DrawingInfo
        DrawingInfo class contains information necessary for drawing/plotting the fischer projection such as
        input_smiles, canvas_width_pixels, canvas_height_pixels, iupac_name etc.

    Returns
    -------
    List[str]
    List of the substituents in the chemical notation for Fischer projection plot.

    """

    return [
        prepare_strings_for_fischer_projection_plot(
            c, drawing_info.canvas_width_pixels, drawing_info.canvas_height_pixels
        )
        for c in substituents_condensed_form
    ]
