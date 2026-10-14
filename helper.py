import numpy as np
from pathlib import Path
from config import MODEL_INFORMATION, ENERGY_GRID, ARGON_CHANNELS

EFFIC = Path.home()/"software"/"snowglobes"/"effic"

def calculate_iron_mass(iron_thickness): #in mm
    iron_density = 7.87 * 1000 #in kg / m^3
    iron_thickness = iron_thickness / 1000 #in m

    argon_height = 14
    argon_width = 15.1
    argon_length = 62

    argon_volume = argon_height * argon_length * argon_width

    outer_height = argon_height + 2 * iron_thickness
    outer_width = argon_width + 2 * iron_thickness
    outer_length = argon_length + 2 * iron_thickness

    outer_volume = outer_height * outer_width * outer_length

    iron_volume = outer_volume - argon_volume
    iron_mass = iron_density * iron_volume #in kg

    return iron_mass / 1000 #in t


def calculate_effective_iron_thickness(mass_t): #in t
    mass_kg = 1000 * mass_t
    iron_density = 7.87 * 1000 #in kg / m^3

    argon_height = 14
    argon_width = 15.1
    argon_length = 62

    iron_volume = mass_kg / iron_density

    cubic_coef = 8
    quadratic_coef = 4*(argon_length + argon_height + argon_width)
    linear_coef = 2*(argon_height*argon_width + argon_length*argon_width+argon_length*argon_height)
    const_coef =  - iron_volume
    coefficients = [cubic_coef, quadratic_coef, linear_coef, const_coef]  

    roots = np.roots(coefficients)

    positive_real_root = [root.real for root in roots if root.imag == 0 and root.real > 0]
    if len(positive_real_root) != 1:
        raise ValueError("Unexpected roots: ", roots)

    return positive_real_root[0] * 1000 # in mm



def apply_snowglobes_threshold(detector, material, threshold):
    
    if material == "iron":
        lookup = MODEL_INFORMATION
    elif material == "argon":
        lookup = ARGON_CHANNELS

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    values = (bin_centers > threshold).astype(float)

    for channel in lookup.keys():

        file = EFFIC / detector / f"effic_{channel}_{detector}.dat"

        with open(file, "w") as f:
            f.write("{" + ",".join(f"{x:.1f}" for x in values) + "}\n")



if __name__ == "__main__":
    #apply_snowglobes_threshold("ar40kt", "argon", 0)
    print(calculate_effective_iron_thickness(2038.1))