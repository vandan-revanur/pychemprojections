from math import ceil
import matplotlib.pyplot as plt
import re
from pychemprojections.utils.rdkit_utils import smiles_to_condensed_form, cleanup_Hs
from pychemprojections.utils.logger_utils import get_module_logger

logger = get_module_logger(__name__)


def split_into_chunks(inp_str, text_width_pixels_coords, threshold_pixels):
    inp_str_len = len(inp_str)
    pixels_per_char = int(text_width_pixels_coords // inp_str_len)
    threshold_chars = int(threshold_pixels // pixels_per_char)

    n_chunks = ceil(inp_str_len / threshold_chars)

    s_idx = 0
    e_idx = inp_str[s_idx : s_idx + threshold_chars].rfind("_") + 2

    temp_list = []
    for i in range(n_chunks):
        temp_list.append(inp_str[s_idx:e_idx] + "$")
        s_idx = e_idx

        if s_idx + threshold_chars >= inp_str_len:
            temp_list.append(inp_str[s_idx:])
            break
        else:
            e_idx = inp_str[: s_idx + threshold_chars].rfind("_") + 2

    ret_string = "-\n-$".join(temp_list)

    # If the splitting results in a string where the last chunk is just -\n-$ then remove that
    if ret_string[-4:] == "-\n-$":
        ret_string = ret_string[:-4]

    return ret_string


def get_width_of_text_string(canvas_width, canvas_height, condensed_string):
    DPI = 300
    fig, ax1 = plt.subplots(figsize=(canvas_width / DPI, canvas_height / DPI))
    r = fig.canvas.get_renderer()
    t = plt.text(0, 0, condensed_string)

    bb = t.get_window_extent(renderer=r)
    text_width_pixels_coords = bb.width
    plt.close(fig)

    return text_width_pixels_coords


def prepare_strings_for_fischer_projection_plot(
    condensed_string, canvas_width, canvas_height
):
    regex = "[0-9]+"  # regex for removing cyclic numbering
    wrap_threshold = 0.2
    threshold_pixels = canvas_width * wrap_threshold
    text_width_pixels_coords = get_width_of_text_string(
        canvas_width, canvas_height, condensed_string
    )

    pos_of_numbers = [
        sub_str.start() for sub_str in re.finditer(regex, condensed_string)
    ]
    condensed_list = list(condensed_string)

    for index in sorted(pos_of_numbers, reverse=True):
        condensed_list.insert(index, "_")

    fischer_projection_notation = "".join(condensed_list)

    if text_width_pixels_coords >= threshold_pixels:
        logger.info(
            f"wrapping text since the length of text is exceeding the {(wrap_threshold + 0.5) * 100}% mark of canvas width"
        )
        fischer_projection_notation = split_into_chunks(
            fischer_projection_notation, text_width_pixels_coords, threshold_pixels
        )

    if fischer_projection_notation[-1] != "$":
        fischer_projection_notation = f"${fischer_projection_notation}$"
    else:
        fischer_projection_notation = f"${fischer_projection_notation}"

    return fischer_projection_notation


def get_condensed_smiles_of_all_substituents(substituents_with_Hs_added):
    return [
        smiles_to_condensed_form(cleanup_Hs(compound))
        for compound in substituents_with_Hs_added
    ]
