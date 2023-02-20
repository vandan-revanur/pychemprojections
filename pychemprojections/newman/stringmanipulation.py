import re
from pychemprojections.utils.rdkit_utils import (
    convert_smiles_to_molecular_formula,
    smiles_to_condensed_form,
    cleanup_Hs,
)
from pychemprojections.utils.smiles_manipulation_utils import (
    get_index_of_corresponding_bracket,
)
from pychemprojections.utils.rdkit_utils import get_atom_ids_of_all_carbon_atoms
from rdkit import Chem


def get_condensed_groups(groups):
    front_and_rear_groups = []
    for k, v in groups.items():
        front_and_rear_groups.extend(v.values())

    groups_condensed = [
        smiles_to_condensed_form(cleanup_Hs(grp)) for grp in front_and_rear_groups
    ]

    return groups_condensed


def prepare_strings_for_newman_projection_plot(condensed_string):
    regex = r"[0-9]+"
    end_pos_of_numbers = sorted(
        [sub_str.end() for sub_str in re.finditer(regex, condensed_string)],
        reverse=True,
    )
    condensed_list = list(condensed_string)

    for e_index in end_pos_of_numbers:
        condensed_list.insert(e_index, "}")

    condensed_string = "".join(condensed_list)
    regex = r"[0-9]+"
    start_pos_of_numbers = sorted(
        [sub_str.start() for sub_str in re.finditer(regex, condensed_string)],
        reverse=True,
    )

    for s_index in start_pos_of_numbers:
        condensed_list.insert(s_index, "_{")

    newman_projection_notation = "".join(condensed_list)
    newman_projection_notation = f"${newman_projection_notation}$"
    return newman_projection_notation


def separate_atoms_and_nums(condensed_string):
    regex = r"[0-9]+"
    end_pos_of_numbers = sorted(
        [sub_str.end() for sub_str in re.finditer(regex, condensed_string)],
        reverse=True,
    )
    condensed_list = list(condensed_string)

    for e_index in end_pos_of_numbers:
        condensed_list.insert(e_index, "_")

    condensed_string = "".join(condensed_list)
    regex = r"[A-Z]+"
    start_pos_of_numbers = sorted(
        [sub_str.end() for sub_str in re.finditer(regex, condensed_string)],
        reverse=True,
    )

    for s_index in start_pos_of_numbers:
        condensed_list.insert(s_index, "_")

    molecular_formula_atoms_nums_separated = "".join(condensed_list)

    return molecular_formula_atoms_nums_separated


def get_c_ids_in_smiles(smiles_mol_prepared):
    carbon_atom_ids_in_smiles = []

    for i in range(len(smiles_mol_prepared) - 1):
        if smiles_mol_prepared[i] == "C" and not smiles_mol_prepared[i + 1].islower():
            carbon_atom_ids_in_smiles.append(i)
    return carbon_atom_ids_in_smiles


def recalibrate_mol_formula_to_account_for_H_addition(molecular_formula):
    atoms_and_nums_info = separate_atoms_and_nums(molecular_formula).split("_")
    num_Hs_idx = atoms_and_nums_info.index("H") + 1
    atoms_and_nums_info[num_Hs_idx] = str(int(atoms_and_nums_info[num_Hs_idx]) - 1)
    molecular_formula = "".join(atoms_and_nums_info)
    return molecular_formula


def smiles_post_processing(input_smiles, condensed, n_carbons_for_truncation):
    if len(get_c_ids_in_smiles(input_smiles)) >= n_carbons_for_truncation:
        molecular_formula = convert_smiles_to_molecular_formula(input_smiles)
        molecular_formula = recalibrate_mol_formula_to_account_for_H_addition(
            molecular_formula
        )
        post_processed_smiles = molecular_formula
    else:
        post_processed_smiles = condensed

    return post_processed_smiles


def get_substituents_front_and_rear(
    smiles_mol_prepared, start_carbon_atom_idx_in_str, end_carbon_atom_idx_in_str
):
    groups = {
        "front": {"1": None, "2": None, "3": None},
        "rear": {"1": None, "2": None, "3": None},
    }
    start_idx_group_1_front = 0
    end_idx_group_1_front = start_carbon_atom_idx_in_str
    group_1_front = smiles_mol_prepared[start_idx_group_1_front:end_idx_group_1_front]

    start_idx_group_2_front = end_idx_group_1_front + 1
    end_idx_group_2_front = get_index_of_corresponding_bracket(
        smiles_mol_prepared,
        start_idx_group_2_front,
        start_bracket_type="(",
        end_bracket_type=")",
    )
    group_2_front = smiles_mol_prepared[
        start_idx_group_2_front + 1 : end_idx_group_2_front
    ]

    start_idx_group_3_front = end_idx_group_2_front + 1
    end_idx_group_3_front = get_index_of_corresponding_bracket(
        smiles_mol_prepared,
        start_idx_group_3_front,
        start_bracket_type="(",
        end_bracket_type=")",
    )
    group_3_front = smiles_mol_prepared[
        start_idx_group_3_front + 1 : end_idx_group_3_front
    ]

    groups["front"]["1"] = group_1_front
    groups["front"]["2"] = group_2_front
    groups["front"]["3"] = group_3_front

    start_idx_group_1_rear = end_carbon_atom_idx_in_str + 1
    end_idx_group_1_rear = get_index_of_corresponding_bracket(
        smiles_mol_prepared,
        start_idx_group_1_rear,
        start_bracket_type="(",
        end_bracket_type=")",
    )
    group_1_rear = smiles_mol_prepared[
        start_idx_group_1_rear + 1 : end_idx_group_1_rear
    ]

    start_idx_group_2_rear = end_idx_group_1_rear + 1
    end_idx_group_2_rear = get_index_of_corresponding_bracket(
        smiles_mol_prepared,
        start_idx_group_2_rear,
        start_bracket_type="(",
        end_bracket_type=")",
    )
    group_2_rear = smiles_mol_prepared[
        start_idx_group_2_rear + 1 : end_idx_group_2_rear
    ]

    start_idx_group_3_rear = end_idx_group_2_rear + 1
    group_3_rear = smiles_mol_prepared[start_idx_group_3_rear:]

    groups["rear"]["1"] = group_1_rear
    groups["rear"]["2"] = group_2_rear
    groups["rear"]["3"] = group_3_rear

    return groups


def get_post_processed_smiles(
    smiles_groups, groups_condensed, n_carbons_for_truncation
):
    post_processed_smiles = []
    for input_smiles, condensed_group in zip(smiles_groups, groups_condensed):
        post_processed_smiles.append(
            smiles_post_processing(
                input_smiles, condensed_group, n_carbons_for_truncation
            )
        )
    return post_processed_smiles


def create_mapping_between_atom_ids_in_smiles_and_rdkit_mol_def(
    mol, smiles_mol_prepared
):
    carbon_atom_ids_in_mol = get_atom_ids_of_all_carbon_atoms(mol)
    carbon_atom_ids_in_str = get_c_ids_in_smiles(smiles_mol_prepared)
    assert len(carbon_atom_ids_in_str) == len(
        carbon_atom_ids_in_mol
    ), "number of carbon atoms in SMILES string not matching the number of carbons in rdkit molecule definition"
    carbon_ids_mol_to_str_map = {
        c_idx_mol: c_idx_str
        for c_idx_str, c_idx_mol in zip(carbon_atom_ids_in_str, carbon_atom_ids_in_mol)
    }

    return carbon_ids_mol_to_str_map


def calibration_for_post_processing(groups):
    smiles_groups = []
    group_1_front_mol = Chem.MolFromSmiles(groups["front"]["1"] + "[H]")
    group_1_front_mol = Chem.RemoveHs(group_1_front_mol)
    smiles_mol_group_1_front = Chem.MolToSmiles(group_1_front_mol)
    smiles_groups.append(smiles_mol_group_1_front)

    for group in [
        groups["front"]["2"],
        groups["front"]["3"],
        groups["rear"]["1"],
        groups["rear"]["2"],
        groups["rear"]["3"],
    ]:
        if "C" in group and "([H])" in group:
            mol_group = Chem.MolFromSmiles("[H]" + group)
            mol_group = Chem.RemoveHs(mol_group)
            smiles_mol_group = Chem.MolToSmiles(mol_group)
        else:
            smiles_mol_group = group
        smiles_groups.append(smiles_mol_group)

    return smiles_groups
