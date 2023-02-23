from pychemprojections.newman.visualization import (
    get_newman_drawing_info,
    plot_newman_projection,
)
import argparse

from pychemprojections.utils.logger_utils import get_module_logger

logger = get_module_logger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--canvas-width-pixels", dest="canvas_width_pixels", required=True, type=int
)
parser.add_argument(
    "--canvas-height-pixels", dest="canvas_height_pixels", required=True, type=int
)
parser.add_argument("--input-smiles", dest="input_smiles", required=True, type=str)
parser.add_argument(
    "--start-carbon-idx", dest="start_carbon_idx", required=True, type=int
)
parser.add_argument("--end-carbon-idx", dest="end_carbon_idx", required=True, type=int)

args = parser.parse_args()
input_smiles = args.input_smiles
canvas_width_pixels = args.canvas_width_pixels
canvas_height_pixels = args.canvas_height_pixels
start_carbon = args.start_carbon_idx
end_carbon = args.end_carbon_idx

logger.info(f"input SMILES: {input_smiles}")
logger.info(f"canvas_width_pixels: {canvas_width_pixels}")
logger.info(f"canvas_height_pixels: {canvas_height_pixels}")

newman_drawing_info = get_newman_drawing_info(
    input_smiles, canvas_width_pixels, canvas_height_pixels, (start_carbon, end_carbon)
)
plot_newman_projection(newman_drawing_info)
