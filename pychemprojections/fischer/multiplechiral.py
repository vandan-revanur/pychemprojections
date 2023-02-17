from pychemprojections.utils.smiles_manipulation_utils import (
    get_index_of_corresponding_bracket,
)
from pychemprojections.fischer.stringmanipulation import (
    get_condensed_smiles_of_all_substituents,
)
from pychemprojections.fischer.singlechiral import (
    get_fisher_notation_of_all_substituents_single_chiral,
)
from typing import List, Dict, Any


def get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
    smiles_mol_prepared: str,
    chiral_tag_of_c_atom: bool = True,
    sq_bracket_begin_first_chiral_center: int = None,
    normal_bracket_begin_achiral_center: int = None,
) -> List[str]:
    if chiral_tag_of_c_atom:
        sq_bracket_begin = sq_bracket_begin_first_chiral_center
        sq_bracket_end = get_index_of_corresponding_bracket(
            smiles_mol_prepared,
            sq_bracket_begin,
            start_bracket_type="[",
            end_bracket_type="]",
        )

        if "@H]" in smiles_mol_prepared:
            sub_group_a = smiles_mol_prepared[:sq_bracket_begin]
            sub_group_b = "H"
            mol_c_normal_bracket_begin_idx = sq_bracket_end + 1
            mol_c_normal_bracket_end_idx = get_index_of_corresponding_bracket(
                smiles_mol_prepared, mol_c_normal_bracket_begin_idx
            )
            sub_group_c = smiles_mol_prepared[
                mol_c_normal_bracket_begin_idx + 1 : mol_c_normal_bracket_end_idx
            ]
            sub_group_d = smiles_mol_prepared[mol_c_normal_bracket_end_idx + 1 :]

        else:
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

        substituent_groups = [sub_group_a, sub_group_b, sub_group_c, sub_group_d]

    else:
        normal_bracket_begin = normal_bracket_begin_achiral_center
        normal_bracket_end = get_index_of_corresponding_bracket(
            smiles_mol_prepared,
            normal_bracket_begin,
            start_bracket_type="(",
            end_bracket_type=")",
        )

        sub_group_a = smiles_mol_prepared[normal_bracket_begin + 1 : normal_bracket_end]

        normal_bracket_begin = normal_bracket_end + 1
        normal_bracket_end = get_index_of_corresponding_bracket(
            smiles_mol_prepared,
            normal_bracket_begin,
            start_bracket_type="(",
            end_bracket_type=")",
        )

        sub_group_b = smiles_mol_prepared[normal_bracket_begin + 1 : normal_bracket_end]

        substituent_groups = [sub_group_a, sub_group_b]

    return substituent_groups


def get_carbon_neighbours_info_multiple_chiral(
    smiles_mol_prepared: str,
) -> Dict[str:Any]:
    smiles_mol_prepared = smiles_mol_prepared.replace("[H]", "H")
    sq_bracket_begin_first_chiral_center = smiles_mol_prepared.find("[")
    sq_bracket_end_last_chiral_center = smiles_mol_prepared.rfind("]")

    carbons_neighbours_info = []
    chiral_config = None
    for i in range(
        sq_bracket_begin_first_chiral_center, sq_bracket_end_last_chiral_center + 1
    ):
        if smiles_mol_prepared[i] == "[":
            subs = get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
                smiles_mol_prepared, sq_bracket_begin_first_chiral_center=i
            )
            carbon_atom_pos_in_string = i + 1
            carbon_type = "chiral"
            if smiles_mol_prepared[i + 2 : i + 4] == "@@":
                chiral_config = "S"
            else:
                chiral_config = "R"

            carbons_neighbours_info.append(
                {
                    "carbon_type": carbon_type,
                    "carbon_atom_pos_in_string": carbon_atom_pos_in_string,
                    "substituents": subs,
                    "config": chiral_config,
                }
            )

    i = sq_bracket_begin_first_chiral_center
    while i <= sq_bracket_end_last_chiral_center:
        if smiles_mol_prepared[i : i + 2] == ")C":
            subs = get_smiles_of_chiral_substituent_groups_in_multiple_chiral_chain(
                smiles_mol_prepared,
                chiral_tag_of_c_atom=False,
                normal_bracket_begin_achiral_center=i + 2,
            )

            carbon_atom_pos_in_string = i + 1
            carbon_type = "achiral"
            carbons_neighbours_info.append(
                {
                    "carbon_type": carbon_type,
                    "carbon_atom_pos_in_string": carbon_atom_pos_in_string,
                    "substituents": subs,
                    "config": chiral_config,
                }
            )
            i = i + 2 + len(subs[0]) + 1 + len(subs[1]) + 1
        else:
            i = i + 1

    return carbons_neighbours_info


def get_right_and_left_substituents(
    substituents_condensed_form: List[Dict[str:Any]],
) -> tuple[List[str], List[str]]:
    right = []
    left = []

    # Arranging the substituents for the first chiral center
    first_group_info = substituents_condensed_form[0]
    substituents = first_group_info["substituents"]
    if first_group_info["config"] == "R":
        right.append(substituents[1])
        left.append(substituents[2])

    else:
        right.append(substituents[2])
        left.append(substituents[1])

    for c_info in substituents_condensed_form[1:-1]:
        substituents = c_info["substituents"]
        if c_info["carbon_type"] == "achiral":
            right.append(substituents[0])
            left.append(substituents[1])

        else:
            if c_info["config"] == "R":
                right.append(substituents[0])
                left.append(substituents[1])
            else:
                right.append(substituents[1])
                left.append(substituents[0])

    # Arranging the substituents for the last chiral center
    last_group_info = substituents_condensed_form[-1]
    substituents = last_group_info["substituents"]
    if last_group_info["config"] == "R":
        right.append(substituents[0])
        left.append(substituents[2])

    else:
        right.append(substituents[2])
        left.append(substituents[0])

    return left, right


def get_condensed_form_info_of_substituents_multiple_chiral(
    sorted_carbons_neighbours_info: List[Dict[str:Any]],
) -> List[Dict[str:Any]]:
    substituents_condensed_form = []

    for c_info in sorted_carbons_neighbours_info:
        c_info["substituents"] = get_condensed_smiles_of_all_substituents(
            c_info["substituents"]
        )
        substituents_condensed_form.append(c_info)

    return substituents_condensed_form


def get_fisher_notation_info_for_substituents_multiple_chiral(
    substituents_condensed_form: List[Dict[str:Any]],
    canvas_width: int,
    canvas_height: int,
) -> List[Dict[str:Any]]:
    substituents_fischer_notation = []

    for c_info in substituents_condensed_form:
        c_info["substituents"] = get_fisher_notation_of_all_substituents_single_chiral(
            c_info["substituents"], canvas_width, canvas_height
        )
        substituents_fischer_notation.append(c_info)

    return substituents_fischer_notation


def remove_unnecessary_neighbours(
    sorted_carbons_neighbours_info: List[Dict[str:Any]],
) -> List[Dict[str:Any]]:
    sorted_carbons_neighbours_info[0]["substituents"].pop(-1)
    sorted_carbons_neighbours_info[-1]["substituents"].pop(0)

    for i in range(0, len(sorted_carbons_neighbours_info)):
        if i != 0 and i != len(sorted_carbons_neighbours_info) - 1:
            if sorted_carbons_neighbours_info[i]["carbon_type"] == "chiral":
                sorted_carbons_neighbours_info[i]["substituents"].pop(0)
                sorted_carbons_neighbours_info[i]["substituents"].pop(-1)

    return sorted_carbons_neighbours_info
