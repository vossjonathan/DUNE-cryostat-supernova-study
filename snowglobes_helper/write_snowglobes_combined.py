import numpy as np
import uproot

from config import DATA_LOOKUP, ENERGY_GRID, ENERGY_BIN_CENTERS, SNOWGLOBES_DIR

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


def build_smearing_matrix(channel, layer):

    bin_centers = ENERGY_BIN_CENTERS

    channel_config = DATA_LOOKUP[layer]["channels"][channel]
    threshold = channel_config["scan_min_energy"]

    smearing_matrix = np.zeros((len(ENERGY_BIN_CENTERS),) * 2)

    for i, E in enumerate(bin_centers):

        if E < threshold:
            continue

        folder = channel_config["folder"]
        filename = channel_config["filename"].format(E=E)

        file = folder/filename

        if not file.exists():
            raise FileNotFoundError(f"{channel}: missing file at E = {E:.5f} MeV:\n{file}")

        smearing_matrix[:, i] = get_smearing_column(file)

    validate_matrix(smearing_matrix)

    return smearing_matrix

def get_detector(layer):

    if layer == "PrimaryIron":
        return "FD_cryostat_primary"
    elif layer == "OuterIron":
        return "FD_cryostat_outer"
    elif layer == "PU":
        return "FD_cryostat_PU"
    else:
        raise ValueError("Not a valid layer")

def write_snowglobes_matrix(smearing_matrix, channel, layer):

    validate_matrix(smearing_matrix)
    n_reco, n_true = smearing_matrix.shape

    detector = get_detector(layer)

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
    matrix = build_smearing_matrix("nue_Fe56", "PrimaryIron")
    write_snowglobes_matrix(matrix, "nue_Fe56", "PrimaryIron")

    # Generate response matrices for all layers and channels
    #for layer in ["PrimaryIron", "PU", "OuterIron"]:
    #    for channel in DATA_LOOKUP[layer]["channels"]:
    #        print(f"Building matrix for {layer}, {channel}")
    #        matrix = build_smearing_matrix(channel, layer)
    #        write_snowglobes_matrix(matrix, channel, layer)