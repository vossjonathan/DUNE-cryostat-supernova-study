import matplotlib.pyplot as plt
import numpy as np
from config import DATA_LOOKUP, ENERGY_GRID
from pathlib import Path
import uproot
from matplotlib.colors import LogNorm

OUT = Path.home()/"software"/"snowglobes"/"out"
ROOT_DATA = Path.home()/"software"/"geant4"/"geant4-work"/"FDSteel"/"build"/"root_data"
SNOWGLOBES = Path.home()/"software"/"snowglobes"


#ANALYSIS FUNCTIONS
def load_events_data(flux, detector, channel, smeared, suffix = ""):
    if smeared:
        file = OUT/f"{flux}_{channel}_{detector}_events_smeared_unweighted{suffix}.dat"
        data = np.loadtxt(file, max_rows = 200)
        events = data[:, 1]
        return events


    else:
        file = OUT/f"{flux}_{channel}_{detector}_events_unweighted.dat"
        data = np.loadtxt(file, max_rows = 200)
        events = data[:, 1]
        print(events.shape)
        return events
    
def get_detector(layer):
    if layer == "PrimaryIron":
        return "FD_cryostat_primary"
    elif layer == "OuterIron":
        return "FD_cryostat_outer"
    elif layer == "PU":
        return "FD_cryostat_PU"
    else:
        raise ValueError("Not a valid layer")
    
def get_root_file(channel, layer, energy):
    channel_info = DATA_LOOKUP[layer]["channels"][channel]

    folder = channel_info["folder"]
    filename = channel_info["filename"].format(E = energy)
    return folder / filename

def load_Edep(channel, layer, energy):
    root_file = get_root_file(channel, layer, energy)

    with uproot.open(root_file) as file:
        Edep = file["FDSteel"]["EdepLAr"].array(library="np")

    return Edep

def get_single_energy_distribution(channel, layer, energy, dep_edges):
    Edep = load_Edep(channel, layer, energy)
    n_generated = len(Edep)

    counts, _ = np.histogram(Edep, bins = dep_edges)
    column = counts  / n_generated
    return column


def build_distribution_matrix(channel, layer):
    channel_info = DATA_LOOKUP[layer]["channels"][channel]
    scan_min_energy = channel_info["scan_min_energy"]

    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    dep_edges = ENERGY_GRID

    n_true_bins = len(energies)
    n_dep_bins = len(dep_edges) - 1
    distribution_matrix = np.zeros((n_dep_bins, n_true_bins))

    print("dep_edges: ", dep_edges)
    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue
        column = get_single_energy_distribution(channel, layer, E, dep_edges)
        distribution_matrix[:, i] = column

    return distribution_matrix, dep_edges

#FDSteel_1.2mm == PrimaryIron mass of 38,1t
#FDSteel_62.7mm == OuterIron mass of 2,000t
def plot_unsmeared_events_iron(flux, detector, layer):
    total_events = 0
    for channel in DATA_LOOKUP[layer]["channels"].keys():
        events_load = load_events_data(flux, detector, channel, smeared = False)
        weight =  DATA_LOOKUP[layer]["channels"][channel]["weight"]
        events_weighted = events_load * weight
        total_events += sum(events_weighted)
        plt.stairs(events_weighted, ENERGY_GRID, label = f"{channel}")

    plt.legend()
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel("Events per 0.5 MeV bin")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"Unsmeared spectrum {detector}, {flux}, {iron_thickness}mm")
    plt.text(x = 0.71, y = 0.5, s = f"Total: {total_events:.2f}", transform = plt.gca().transAxes, ha = "left", va = "top")
    plt.savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


def plot_distribution_matrix(channel, layer):

    matrix, dep_edges = build_distribution_matrix(channel, layer)
    print(matrix.sum(axis = 0))
    plt.imshow(matrix, origin = "lower", aspect = "auto",    extent=[ENERGY_GRID[0],ENERGY_GRID[-1],dep_edges[0],dep_edges[-1]], norm=LogNorm(vmin=1e-5, vmax=matrix.max()))
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel(r"LAr deposited energy [MeV]")
    plt.title(f"LAr energy deposition distribution for {channel}, {layer}")
    plt.colorbar(label = "Probability")
    plt.savefig(f"plotting/energy_distribution_matrix_{channel}_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_hist_neutrino_energy(channel, layer, energy):

    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    energy = round(energies[np.argmin(np.abs(energies - energy))], 5)

    Edep = load_Edep(channel, layer, energy)

    Ar_dep_energy_grid = np.arange(0, 100.5, 0.5)

    zero_fraction = np.mean(Edep == 0)
    positive_fraction = np.mean(Edep > 0)

    print(f"E_nu = {energy}, exact zero:", np.sum(Edep == 0))
    print(f"E_nu = {energy}, 0 < Edep < 0.5:", np.sum((Edep > 0) & (Edep < 0.5)))
    print(f"E_nu = {energy} Edep >= 0.5:", np.sum(Edep >= 0.5))

    print(f"Exactly zero: {zero_fraction:.4f}")
    print(f"Edep > 0:     {positive_fraction:.4f}")
    plt.hist(Edep, bins = Ar_dep_energy_grid)
    plt.title(fr"$E_{{dep}}$ in LAr distribution for $E_\nu$ = {energy} MeV, {channel}, {layer}")
    plt.xlabel("Edep [MeV]")
    plt.ylabel("Event counts")
    plt.savefig(f"plotting/{energy}_MeV_distribution_{channel}_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def calculate_efficiencies(channel, layer, dep_threshold):
    scan_min_energy = DATA_LOOKUP[layer]["channels"][channel]["scan_min_energy"]

    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    efficiencies = np.zeros(len(energies))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue

        Edep = load_Edep(channel, layer, E)

        efficiencies[i] = np.mean(Edep > dep_threshold)

    return efficiencies


def plot_deposition_efficiency(channel, layer, dep_threshold):

    scan_min_energy = DATA_LOOKUP[layer]["channels"][channel]["scan_min_energy"]
    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    mask = energies >= scan_min_energy
    efficiencies = calculate_efficiencies(channel, layer, dep_threshold)

    plt.plot(energies[mask], efficiencies[mask])
    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.grid(True, alpha = 0.3)
    plt.title(f"Event fraction with deposition energy > {dep_threshold} MeV, {channel}, {layer}")
    plt.savefig(f"plotting/effiency_fraction_{channel}_{layer}_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


def get_layer(detector):
    if detector.split("_")[-1] == "primary":
        return"PrimaryIron"

    elif detector.split("_")[-1] == "outer":
        return"OuterIron"

    elif detector.split("_")[-1] == "PU":
        return"PU"

    else:
        raise ValueError("Fehler bei Detektoreingabe, das sind die validen Detektornamen: FD_cryostat_primary, FD_cryostat_outer, FD_cryostat_PU")



def plot_unsmeared_events_iron(flux, detector):
    total_events = 0
    layer = get_layer(detector)
    channel_config = DATA_LOOKUP[layer]["channels"]
    for channel in channel_config.keys():
        events_load = load_events_data(flux, detector, channel, smeared = False)
        weight = channel_config[channel]["weight"]
        events_weighted = events_load * weight
        total_events += sum(events_weighted)
        plt.stairs(events_weighted, ENERGY_GRID, label = f"{channel}")

    plt.legend()
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel("Events per 0.5 MeV bin")
    plt.grid(True, alpha = 0.3)
    plt.title(f"Unsmeared spectrum {detector}, {flux}, {layer}")
    plt.text(x = 0.71, y = 0.5, s = f"Total: {total_events:.2f}", transform = plt.gca().transAxes, ha = "left", va = "top")
    plt.savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()



def plot_one_channel_weighted(flux, detector, channel, smeared):
    layer = get_layer(detector)

    channel_config = DATA_LOOKUP[layer]["channels"]

    events_load = load_events_data(flux, detector, channel, smeared)
    weight = channel_config[channel]["weight"]
    events_weighted = events_load * weight

    word = "smeared" if smeared else "unsmeared"
    color = "blue" if smeared else "red"

    plt.stairs(events_weighted, ENERGY_GRID, label = word + f" {channel}", color = color)
    plt.text(0.6, 0.86,f"Events {word}: {events_weighted.sum():.4f}",transform=plt.gca().transAxes)

    plt.xlim(0, 50)
    plt.xlabel(r"E [MeV]")
    plt.ylabel("Number of events per bin")
    plt.legend()
    plt.grid(True, alpha = 0.3)
    plt.title(f"{word.capitalize()} events for {detector}, {channel}, {flux} flux")
    plt.savefig(f"plotting/{flux}_{detector}_{channel}_{word}_events.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_one_channel_smeared_vs_unsmeared(flux, channel, layer):

    detector = get_detector(layer)

    channel_config = DATA_LOOKUP[layer]["channels"][channel]

    weight = channel_config["weight"]

    events_unsmeared = load_events_data(flux, detector, channel, smeared = False) * weight

    events_smeared = load_events_data(flux, detector, channel, smeared = True) * weight

    threshold = 5
    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    events_smeared_threshold = events_smeared[bin_centers > threshold]


    plt.stairs(events_smeared, ENERGY_GRID, label = "Smeared: Reconstructed E", color = "blue")
    plt.stairs(events_unsmeared, ENERGY_GRID, label = "Unsmeared: True E", color = "red")

    plt.text(0.52, 0.71,f"Total events unsmeared: {events_unsmeared.sum():.2f}\nTotal events smeared: {events_smeared.sum():.2f}\nEvents with > 5 MeV: {events_smeared_threshold.sum():.2f}",transform=plt.gca().transAxes)

    plt.xlabel(r" E [MeV]")
    plt.xlim(0, 50)
    plt.ylabel("Number of Events per bin")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.title(f"{detector}, {channel} unsmeared vs. smeared, {flux} flux")
    plt.savefig(f"plotting/{flux}_{channel}_{detector}_smeared_vs_unsmeared_new.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_all_iron_smeared_vs_unsmeared(flux, layer):
    channel_config = DATA_LOOKUP[layer]["channels"]
    detector = get_detector(layer)
    eES_events_smeared = np.zeros(200)
    eES_events_unsmeared = np.zeros(200)
    for channel in channel_config.keys():
        if channel.endswith("e"):
            weight = channel_config[channel]["weight"]
            events_smeared = load_events_data(flux, detector, channel, smeared = True) * weight
            events_unsmeared = load_events_data(flux, detector, channel, smeared = False) * weight
            eES_events_smeared += events_smeared
            eES_events_unsmeared += events_unsmeared
        else:
            continue

    weight_nue = channel_config["nue_Fe56"]["weight"]
    events_smeared_nue_Fe56 = load_events_data(flux, detector, "nue_Fe56", smeared = True) * weight_nue
    events_unsmeared_nue_Fe56 = load_events_data(flux, detector, "nue_Fe56", smeared = False) * weight_nue
    

    weight_nuebar = channel_config["nuebar_Fe56"]["weight"]
    events_smeared_nuebar_Fe56 = load_events_data(flux, detector, "nuebar_Fe56", smeared = True) * weight_nuebar
    events_unsmeared_nuebar_Fe56 = load_events_data(flux, detector, "nuebar_Fe56", smeared = False) * weight_nuebar

    plt.plot([], [], color="gray", linestyle="--", label="unsmeared events")

    plt.stairs(eES_events_unsmeared, ENERGY_GRID, color = "tab:green", linestyle = "--")
    plt.stairs(eES_events_smeared, ENERGY_GRID, label = "eES", color = "tab:green", linestyle = "-", baseline = None)

    plt.stairs(events_unsmeared_nue_Fe56, ENERGY_GRID, color = "tab:blue", linestyle = "--")
    plt.stairs(events_smeared_nue_Fe56, ENERGY_GRID, label = r"CC $\nu_e$", color = "tab:blue", linestyle = "-", baseline = None)

    plt.stairs(events_unsmeared_nuebar_Fe56, ENERGY_GRID, color = "tab:orange", linestyle = "--")
    plt.stairs(events_smeared_nuebar_Fe56, ENERGY_GRID, label = r"CC $\bar\nu_e$", color = "tab:orange", linestyle = "-", baseline = None)


    total_events_unsmeared = eES_events_unsmeared + events_unsmeared_nue_Fe56 + events_unsmeared_nuebar_Fe56
    total_events_smeared = eES_events_smeared + events_smeared_nue_Fe56 + events_smeared_nuebar_Fe56

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    total_threshold = total_events_smeared[bin_centers > 5]

    plt.stairs(total_events_smeared, ENERGY_GRID, color = "black", linestyle = "-", label = "Total smeared", linewidth = 1.6, baseline = None)

    plt.text(0.56, 0.53,f"Total events unsmeared: {total_events_unsmeared.sum():.2f}\nTotal events smeared: {total_events_smeared.sum():.2f}\nEvents with > 5 MeV: {total_threshold.sum():.2f}",transform=plt.gca().transAxes, linespacing = 1.4)

    plt.xlabel(r"Smeared: E in LAr, Unsmeared: $E_\nu$ [MeV]")
    plt.xlim(0, 50)
    plt.ylabel("Number of Events per bin")
    plt.legend()
    plt.grid(True, alpha=0.3)

    #plt.title(f"All implemented Fe channels, {flux}, {detector}")
    plt.savefig(f"plotting/{flux}__{detector}_all_iron.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


if __name__ == "__main__":
    '''   layers = ["PrimaryIron", "OuterIron"]
    for layer in layers:
        channels = DATA_LOOKUP[layer]["channels"]
        for channel in channels.keys():
            #plot_distribution_matrix(channel, layer)
            plot_deposition_efficiency(channel, layer, 0.0)
            plot_deposition_efficiency(channel, layer, 0.5)
            plot_deposition_efficiency(channel, layer, 5.0)
    '''
    #plot_distribution_matrix("ibd", "PU")
    #plot_deposition_efficiency("ibd", "PU", 0.0)
    #plot_deposition_efficiency("ibd", "PU", 0.4)

    #plot_unsmeared_events_iron("garching", "FD_cryostat_primary")
    #plot_unsmeared_events_iron("garching", "FD_cryostat_outer")
    #plot_one_channel_weighted("garching", "FD_cryostat_PU", "ibd", smeared = False)
    #plot_one_channel_smeared_vs_unsmeared("garching", "ibd", "PU")
    '''    unsmeared = load_events_data(
    "garching", "FD_cryostat_PU", "ibd", False
    )

    for threshold in [0.0, 0.5, 1.0, 3.0, 5.0, 10.0]:

        efficiencies = calculate_efficiencies(
            "ibd", "PU", threshold
        )

        events = np.sum(unsmeared * efficiencies)

        print(
            f"Edep > {threshold:.1f} MeV: "
            f"{events:.4f} events"
        )'''

    #plot_one_channel_smeared_vs_unsmeared("garching", "ibd", "PU")
    '''    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    E = energies[np.argmin(np.abs(energies - 15.0))]

    Edep = load_Edep("ibd", "PU", E)

    bins = np.arange(0.4, 1.1, 0.005)

    plt.hist(Edep, bins=bins)
    plt.axvline(0.511, color="red", linestyle="--", label="511 keV")
    plt.xlabel("LAr deposited energy [MeV]")
    plt.ylabel("Events per 5 keV bin")
    plt.legend()
    plt.show()

    first_bin = (Edep >= 0.5) & (Edep < 0.9975)
    annihilation_region = (Edep >= 0.509) & (Edep < 0.513)

    print("Fraction in first bin:", np.mean(first_bin))
    print(
        "Fraction of first-bin events near 511 keV:",
        np.sum(annihilation_region) / np.sum(first_bin)
    )

    unsmeared = load_events_data(
    "garching", "FD_cryostat_PU", "ibd", False
    )

    smeared = load_events_data(
        "garching", "FD_cryostat_PU", "ibd", True
    )

    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    probabilities = np.zeros(len(energies))

    for i, E in enumerate(energies):
        if E < 1.804:
            continue

        Edep = load_Edep("ibd", "PU", E)

        probabilities[i] = np.mean(
            (Edep >= ENERGY_GRID[0]) &
            (Edep < ENERGY_GRID[1])
        )

    first_bin_direct = np.dot(unsmeared, probabilities)

    print("Direct ROOT calculation:", first_bin_direct)
    print("SNOwGLoBES first bin:", smeared[0])
    print("Difference:", first_bin_direct - smeared[0])'''
    #plot_all_iron_smeared_vs_unsmeared("garching", "PrimaryIron")
    plot_all_iron_smeared_vs_unsmeared("garching", "OuterIron")