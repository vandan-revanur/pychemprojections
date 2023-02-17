from pychemprojections.fischer.visualization import plot_fischer_projection

canvas_width_pixels = 3000
canvas_height_pixels = 3000
input_smiles_for_fischer = "C1CC1[C@](Br)(I)(Cl)"

plot_fischer_projection(
    input_smiles_for_fischer, canvas_width_pixels, canvas_height_pixels
)
