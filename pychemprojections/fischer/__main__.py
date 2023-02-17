from pychemprojections.fischer.visualization import plot_fischer_projection
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

args = parser.parse_args()
input_smiles = args.input_smiles
canvas_width_pixels = args.canvas_width_pixels
canvas_height_pixels = args.canvas_height_pixels

logger.info(f"input SMILES: {input_smiles}")
logger.info(f"canvas_width_pixels: {canvas_width_pixels}")
logger.info(f"canvas_height_pixels: {canvas_height_pixels}")

plot_fischer_projection(input_smiles, canvas_width_pixels, canvas_height_pixels)
