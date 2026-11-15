import numpy as np
from config import MODEL_INFORMATION, ARGON_CHANNELS, SNOWGLOBES_DIR, ENERGY_BIN_CENTERS

EFFIC = SNOWGLOBES_DIR / "effic"

def calculate_iron_mass(iron_thickness): #in mm

    if not np.isfinite(iron_thickness) or iron_thickness <= 0:
        raise ValueError("iron_thickness must be finite and positive")
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

#calculate iron thicknes for the only iron cryostat model
#The iron is directly next to the LAr
def calculate_effective_iron_thickness_onlyiron(mass_t): #in t

    if not np.isfinite(mass_t) or mass_t <= 0:
        raise ValueError("mass_t must be finite and positive")
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

    positive_real_root = [root.real for root in roots if np.isclose(root.imag, 0) and root.real > 0]
    if len(positive_real_root) != 1:
        raise ValueError(f"Expected one positive real physical root; got {roots}")

    return positive_real_root[0] * 1000 # in mm


#calculates thickness of outer iron layer
#uses the more complex cryostat model including a 80cm polyurethane layer between the primary iron membrane and outer iron structure
def calculate_outer_iron_thickness_with_PU(mass_t): #in t

    if not np.isfinite(mass_t) or mass_t <= 0:
        raise ValueError("mass_t must be finite and positive")
    mass_kg = 1000 * mass_t
    iron_density = 7.87 * 1000 #in kg / m^3

    argon_height = 14
    argon_width = 15.1
    argon_length = 62

    primary_iron_thickness = 1.2e-3
    PU_thickness = 0.8

    combined_height = argon_height + 2 * (primary_iron_thickness + PU_thickness)
    combined_width = argon_width + 2 * (primary_iron_thickness + PU_thickness)
    combined_length = argon_length + 2 * (primary_iron_thickness + PU_thickness)

    iron_volume = mass_kg / iron_density

    cubic_coef = 8
    quadratic_coef = 4*(combined_length + combined_height + combined_width)
    linear_coef = 2*(combined_height*combined_width + combined_length*combined_width+combined_length*combined_height)
    const_coef =  - iron_volume
    coefficients = [cubic_coef, quadratic_coef, linear_coef, const_coef]

    roots = np.roots(coefficients)

    positive_real_root = [root.real for root in roots if np.isclose(root.imag, 0) and root.real > 0]
    if len(positive_real_root) != 1:
        raise ValueError(f"Expected one positive real physical root; got {roots}")

    return positive_real_root[0] * 1000 # in mm


def apply_snowglobes_threshold(detector, material, threshold):

    if material == "iron":
        lookup = MODEL_INFORMATION
    elif material == "argon":
        lookup = ARGON_CHANNELS
    else:
        raise ValueError(f"Unsupported material {material!r}; choose iron or argon")
    bin_centers = ENERGY_BIN_CENTERS

    values = (bin_centers > threshold).astype(float)

    for channel in lookup.keys():

        file = EFFIC / detector / f"effic_{channel}_{detector}.dat"

        file.parent.mkdir(parents=True, exist_ok=True)
        with open(file, "w") as f:
            f.write("{" + ",".join(f"{x:.1f}" for x in values) + "}\n")


if __name__ == "__main__":
    print(calculate_outer_iron_thickness_with_PU(2000))
