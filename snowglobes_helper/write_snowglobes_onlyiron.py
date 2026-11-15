import numpy as np
import uproot

from config import MODEL_INFORMATION, ENERGY_GRID, ENERGY_BIN_CENTERS, SNOWGLOBES_DIR

SMEAR = SNOWGLOBES_DIR / "smear"

def get_smearing_column(file):

    with uproot.open(file) as root_file:
        try:
            Edep = root_file["FDSteel"]["EdepLAr"].array(library="np")
        except KeyError as exc:
            raise ValueError(f"{file}: missing FDSteel/EdepLAr") from exc
    n_generated = len(Edep)
    if not n_generated:
        raise ValueError(f"{file}: zero generated events")
    if Edep.ndim != 1 or not np.isfinite(Edep).all():
        raise ValueError(f"{file}: deposited energies must be finite scalar MeV values")
    counts, _ = np.histogram(Edep, bins=ENERGY_GRID)
    return counts / n_generated


def validate_matrix(matrix):

    expected = (len(ENERGY_GRID) - 1,) * 2
    if matrix.shape != expected:
        raise ValueError(f"Expected matrix shape {expected}, got {matrix.shape}")
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError("Matrix must contain finite nonnegative probabilities")
    if (matrix.sum(axis=0) > 1 + 1e-12).any():
        raise ValueError("Matrix probability column sum exceeds one")


def build_smearing_matrix(channel, iron_thickness):

    bin_centers = ENERGY_BIN_CENTERS

    channel_config = MODEL_INFORMATION[channel]
    threshold = channel_config["scan_min_energy"]
    simulation_config = channel_config["simulations"][iron_thickness]

    smearing_matrix = np.zeros((len(ENERGY_BIN_CENTERS),) * 2)

    for i, E in enumerate(bin_centers):

        if E < threshold:
            continue

        folder = simulation_config["folder"]
        filename = simulation_config["filename"].format(E=E)

        file = folder / filename

        if not file.exists():
            raise FileNotFoundError(f"{channel}: missing file at E = {E:.5f} MeV:\n{file}")

        smearing_matrix[:, i] = get_smearing_column(file)

    validate_matrix(smearing_matrix)

    return smearing_matrix


def write_snowglobes_matrix(smearing_matrix, channel, iron_thickness):

    validate_matrix(smearing_matrix)
    n_reco, n_true = smearing_matrix.shape

    detector = f"FDSteel_{iron_thickness}mm"

    filename = SMEAR/f"{detector}"/f"smear_{channel}_{detector}.dat"
    filename.parent.mkdir(parents=True, exist_ok=True)
    with open(filename, "w") as f:

        f.write(f"energy(#{channel}_smear)<\n")

        for i in range(n_reco):

            values = ",".join("0.0" if x == 0.0 else f"{x:.8e}" for x in smearing_matrix[i])
            line = f"{{0,{n_true - 1},{values}}}"

            if i == 0:
                f.write("@energy = " + line + ":\n")

            elif i == n_reco - 1:
                f.write(line + ";\n")

            else:
                f.write(line + ":\n")

        f.write(">\n")

if __name__ == "__main__":

    # Generate response matrix for a single channel
    channel = "nue_Fe56"
    iron_thickness = 12

    matrix = build_smearing_matrix(channel, iron_thickness)
    write_snowglobes_matrix(matrix, channel, iron_thickness)

    # Generate response matrices for all channels and iron thicknesses
    #for iron_thickness in [1.2, 12, 62.7]:
    #    for channel in MODEL_INFORMATION:
    #        print(f"Building matrix for {channel}, {iron_thickness} mm")
    #        matrix = build_smearing_matrix(channel, iron_thickness)
    #        write_snowglobes_matrix(matrix, channel, iron_thickness)