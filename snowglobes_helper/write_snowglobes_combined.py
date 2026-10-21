import numpy as np
import uproot
from pathlib import Path
from config import MODEL_INFORMATION, ENERGY_GRID, DATA_LOOKUP


SNOWGLOBES = Path.home()/"software"/"snowglobes"
SMEAR = SNOWGLOBES/"smear"
ROOT_DATA = Path.home()/"software"/"geant4"/"geant4-work"/"FDSteel"/"build"/"root_data"


def get_smearing_column(file):

    file = uproot.open(file)

    tree = file["FDSteel"]
    Edep= tree["EdepLAr"].array(library = "np")
    n_generated = len(Edep)


    bin_edges = np.linspace(0.5, 100, 201)


    counts, edges = np.histogram(Edep, bins = bin_edges)

    column = counts / n_generated


    return column

def build_smearing_matrix(channel, layer):

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    channel_config = DATA_LOOKUP[layer]["channels"][channel]
    threshold = channel_config["scan_min_energy"]

    smearing_matrix = np.zeros((200, 200))

    for i, E in enumerate(bin_centers):

        if E < threshold:
            continue

        folder = channel_config["folder"]
        filename = channel_config["filename"].format(E=E)
        
        file = folder/filename

        if not file.exists():
            raise FileNotFoundError(f"{channel}: missing file at E = {E:.5f} MeV:\n{file}")

        smearing_matrix[:, i] = get_smearing_column(file)

    #if channel_config["ES"]:
    #    flavor = channel.split("_")[0]
    #    smearing_matrix = create_ES_smearing_matrix(smearing_matrix, flavor)

    if np.isnan(smearing_matrix).any():
        raise ValueError(f"{channel}: matrix contains NaNs")

    if np.any(smearing_matrix < 0):
        raise ValueError(f"{channel}: matrix contains negative values")

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
    n_reco, n_true = smearing_matrix.shape

    detector = get_detector(layer)
    
    filename = SMEAR/f"{detector}"/f"smear_{channel}_{detector}.dat"
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
    # matrix = build_smearing_matrix("ibd", "PU")
    # print("Fraction of events which deposit any energy which land in first bin: ", matrix[0, :] / matrix.sum(axis = 0))
    # print("first bin row: ", matrix[0, :] )
    # print("axis = 0 sum: ", matrix.sum(axis = 0))

    '''    layers = ["PrimaryIron", "PU", "OuterIron"]
    for layer in layers:
        print("Building for: ", layer)
        for channel in DATA_LOOKUP[layer]["channels"].keys():
            print("Building for: ", channel)
            matrix = build_smearing_matrix(channel, layer)
            write_snowglobes_matrix(matrix, channel, layer)'''
