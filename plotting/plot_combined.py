import matplotlib.pyplot as plt
import numpy as np
from config import DATA_LOOKUP, ENERGY_GRID, ARGON_CHANNELS
from pathlib import Path
import uproot
from matplotlib.colors import LogNorm

from config import SNOWGLOBES_DIR, ROOT_DATA, PLOT_DIR, ENERGY_BIN_CENTERS

OUT = SNOWGLOBES_DIR / "out"


def _savefig(filename, **kwargs):

    destination = PLOT_DIR / Path(filename).name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f"Figure already exists: {destination}; choose a new PLOT_DIR")
    plt.savefig(destination, **kwargs)


#ANALYSIS FUNCTIONS
def load_events_data(flux, detector, channel, smeared, suffix = ""):

    ending = f"smeared_unweighted{suffix}" if smeared else "unweighted"
    file = OUT / f"{flux}_{channel}_{detector}_events_{ending}.dat"
    data = np.loadtxt(file, max_rows=len(ENERGY_BIN_CENTERS), ndmin=2)
    if data.shape[0] != len(ENERGY_BIN_CENTERS) or data.shape[1] < 2:
        raise ValueError(f"{file}: expected {len(ENERGY_BIN_CENTERS)} rows and at least two columns")
    if not np.isfinite(data).all() or (data[:, 1] < 0).any():
        raise ValueError(f"{file}: expected finite nonnegative event counts")
    return data[:, 1]


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

    if Edep.ndim != 1 or not len(Edep) or not np.isfinite(Edep).all():
        raise ValueError(f"{root_file}: expected nonempty finite scalar deposited energies in MeV")
    return Edep

def get_single_energy_distribution(channel, layer, energy, dep_edges):

    Edep = load_Edep(channel, layer, energy)
    n_generated = len(Edep)
    if not n_generated or not np.isfinite(Edep).all():
        raise ValueError("Expected nonempty, finite deposited energies in MeV")

    counts, _ = np.histogram(Edep, bins = dep_edges)
    column = counts  / n_generated
    return column


def build_distribution_matrix(channel, layer):

    channel_info = DATA_LOOKUP[layer]["channels"][channel]
    scan_min_energy = channel_info["scan_min_energy"]

    energies = ENERGY_BIN_CENTERS
    dep_edges = ENERGY_GRID

    n_true_bins = len(energies)
    n_dep_bins = len(dep_edges) - 1
    distribution_matrix = np.zeros((n_dep_bins, n_true_bins))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue
        column = get_single_energy_distribution(channel, layer, E, dep_edges)
        distribution_matrix[:, i] = column

    return distribution_matrix, dep_edges

#FDSteel_1.2mm == PrimaryIron mass of 38,1t
#FDSteel_62.7mm == OuterIron mass of 2,000t
def plot_unsmeared_events_iron_for_layer(flux, detector, layer):

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
    _savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


def plot_distribution_matrix(channel, layer, title):

    matrix, dep_edges = build_distribution_matrix(channel, layer)
    print(matrix.sum(axis = 0))
    plt.imshow(matrix, origin = "lower", aspect = "auto",    extent=[ENERGY_GRID[0],ENERGY_GRID[-1],dep_edges[0],dep_edges[-1]], norm=LogNorm(vmin=1e-5, vmax=matrix.max()))
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]", fontsize = 12)
    plt.ylabel(r"LAr deposited energy [MeV]", fontsize = 12)
    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    if title:
        plt.title(f"LAr energy deposition distribution for {channel}, {layer}")
    plt.colorbar().set_label("Probability", fontsize=12)
    _savefig(f"plotting/energy_distribution_matrix_{channel}_{layer}.pdf", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_hist_neutrino_energy(channel, layer, energy):

    energies = ENERGY_BIN_CENTERS
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
    _savefig(f"plotting/{energy}_MeV_distribution_{channel}_{layer}.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def calculate_efficiencies(channel, layer, dep_threshold):

    scan_min_energy = DATA_LOOKUP[layer]["channels"][channel]["scan_min_energy"]

    energies = ENERGY_BIN_CENTERS

    efficiencies = np.zeros(len(energies))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue

        Edep = load_Edep(channel, layer, E)

        efficiencies[i] = np.mean(Edep > dep_threshold)

    return efficiencies


def plot_deposition_efficiency(channel, layer, dep_threshold):

    scan_min_energy = DATA_LOOKUP[layer]["channels"][channel]["scan_min_energy"]
    energies = ENERGY_BIN_CENTERS
    mask = energies >= scan_min_energy
    efficiencies = calculate_efficiencies(channel, layer, dep_threshold)

    plt.plot(energies[mask], efficiencies[mask])
    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.grid(True, alpha = 0.3)
    plt.title(f"Event fraction with deposition energy > {dep_threshold} MeV, {channel}, {layer}")
    _savefig(f"plotting/effiency_fraction_{channel}_{layer}_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
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
    _savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{layer}.svg", bbox_inches = "tight", dpi = 900)
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
    _savefig(f"plotting/{flux}_{detector}_{channel}_{word}_events.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_one_channel_smeared_vs_unsmeared(flux, channel, layer):

    detector = get_detector(layer)

    channel_config = DATA_LOOKUP[layer]["channels"][channel]

    weight = channel_config["weight"]

    events_unsmeared = load_events_data(flux, detector, channel, smeared = False) * weight

    events_smeared = load_events_data(flux, detector, channel, smeared = True) * weight

    threshold = 5
    bin_centers = ENERGY_BIN_CENTERS
    events_smeared_threshold = events_smeared[bin_centers > threshold]


    plt.stairs(events_smeared, ENERGY_GRID, label = "LAr deposited energy", color = "tab:blue")
    plt.stairs(events_unsmeared, ENERGY_GRID, label = r"Unsmeared: True $E_\nu$", color = "tab:red")

    plt.text(0.46, 0.63,f"Total events unsmeared: {events_unsmeared.sum():.2f}\nTotal events smeared: {events_smeared.sum():.2f}\nEvents with > 5 MeV: {events_smeared_threshold.sum():.2f}",transform=plt.gca().transAxes, linespacing = 1.4, fontsize = 12)

    plt.xlabel(r" Energy [MeV]", fontsize = 12)
    plt.xlim(0, 50)
    plt.ylabel("Number of Events per bin", fontsize = 12)
    plt.legend(fontsize = 12)
    plt.grid(True, alpha=0.3)

    #plt.title(f"{detector}, {channel} unsmeared vs. smeared, {flux} flux")
    _savefig(f"plotting/{flux}_{channel}_{detector}_smeared_vs_unsmeared_new.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_all_iron_smeared_vs_unsmeared(flux, layer):

    channel_config = DATA_LOOKUP[layer]["channels"]
    detector = get_detector(layer)
    eES_events_smeared = np.zeros(len(ENERGY_BIN_CENTERS))
    eES_events_unsmeared = np.zeros(len(ENERGY_BIN_CENTERS))
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

    #plt.plot([], [], color="gray", linestyle="--", label="unsmeared events")

    #plt.stairs(eES_events_unsmeared, ENERGY_GRID, color = "tab:green", linestyle = "--")
    plt.stairs(eES_events_smeared, ENERGY_GRID, label = "eES", color = "tab:green", linestyle = "-", baseline = None)

    #plt.stairs(events_unsmeared_nue_Fe56, ENERGY_GRID, color = "tab:blue", linestyle = "--")
    plt.stairs(events_smeared_nue_Fe56, ENERGY_GRID, label = r"$\nu_e$ CC on 56Fe", color = "tab:blue", linestyle = "-", baseline = None)

    #plt.stairs(events_unsmeared_nuebar_Fe56, ENERGY_GRID, color = "tab:orange", linestyle = "--")
    plt.stairs(events_smeared_nuebar_Fe56, ENERGY_GRID, label = r"$\bar\nu_e$ CC on 56Fe", color = "tab:orange", linestyle = "-", baseline = None)


    total_events_unsmeared = eES_events_unsmeared + events_unsmeared_nue_Fe56 + events_unsmeared_nuebar_Fe56
    total_events_smeared = eES_events_smeared + events_smeared_nue_Fe56 + events_smeared_nuebar_Fe56
    bin_centers = ENERGY_BIN_CENTERS

    total_threshold = total_events_smeared[bin_centers > 5]

    plt.stairs(total_events_smeared, ENERGY_GRID, color = "black", linestyle = "-", label = "Total smeared", linewidth = 1.6, baseline = None)

    plt.text(0.47, 0.50,f"Total events unsmeared: {total_events_unsmeared.sum():.2f}\nTotal events smeared: {total_events_smeared.sum():.2f}\nEvents with > 5 MeV: {total_threshold.sum():.2f}",transform=plt.gca().transAxes, linespacing = 1.4, fontsize = 12)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.ylim(bottom = 0)
    plt.axvline(5,color="gray",linestyle="--",linewidth=1.2,alpha=0.8)
    plt.xlabel(r"Energy deposited in LAr [MeV]", fontsize = 12)
    plt.xlim(0, 20)
    plt.ylabel("Number of Events per bin", fontsize = 12)
    plt.legend(fontsize = 12)
    plt.grid(True, alpha=0.3)

    #plt.title(f"All implemented Fe channels, {flux}, {detector}")
    _savefig(f"plotting/{flux}__{detector}_all_iron.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_all_argon_channels_smeared_and_cryostat(flux, detector):

    total_cryostat_events = np.zeros(len(ENERGY_BIN_CENTERS))

    total_PrimaryIron_events = np.zeros(len(ENERGY_BIN_CENTERS))
    total_PU_events = np.zeros(len(ENERGY_BIN_CENTERS))
    total_OuterIron_events = np.zeros(len(ENERGY_BIN_CENTERS))

    layers = ["PrimaryIron", "PU", "OuterIron"]
    for layer in layers:
        channel_config = DATA_LOOKUP[layer]["channels"]
        detector_cryostat = get_detector(layer)
        for channel in channel_config.keys():
            events_smeared = load_events_data(flux, detector_cryostat, channel, True) * channel_config[channel]["weight"]
            if layer == "PrimaryIron":
                total_PrimaryIron_events += 4 * events_smeared
            elif layer == "PU":
                total_PU_events += 4 * events_smeared
            elif layer == "OuterIron":
                total_OuterIron_events += 4 * events_smeared
            total_cryostat_events += 4 * events_smeared

    total_events = np.zeros(len(ENERGY_BIN_CENTERS))
    total_eES_events = np.zeros(len(ENERGY_BIN_CENTERS))

    for channel in ARGON_CHANNELS.keys():
        events_smeared = load_events_data(flux, detector, channel, smeared = True)
        weight = ARGON_CHANNELS[channel]
        events_smeared_w = events_smeared *weight
        if channel.endswith("e"):
            total_eES_events += events_smeared_w
        if channel == "nue_Ar40":
            plt.stairs(events_smeared_w, ENERGY_GRID, label = r"$\nu_e$ CC on LAr", baseline= None)
        if channel == "nuebar_Ar40":
            plt.stairs(events_smeared_w, ENERGY_GRID, label = r"$\bar\nu_e$ CC on LAr", baseline= None)
        total_events += events_smeared_w


    print("unmasked: ", total_events)
    bin_centers = ENERGY_BIN_CENTERS
    mask = bin_centers > 5
    total_events_mask = total_events[mask]
    total_cryostat_mask = total_cryostat_events[mask]
    print("masked: ", total_events_mask)
    print("Cryostat total:", total_cryostat_events.sum())
    print("Cryostat > 5 MeV:", total_cryostat_events[mask].sum())
    print("LAr > 5 MeV:", total_events[mask].sum())

    print(
        "Relative contribution:",
        100 * total_cryostat_events[mask].sum()
        / total_events[mask].sum(),
        "%"
)


    plt.stairs(total_eES_events, ENERGY_GRID, label = "Total eES", baseline= None)

    plt.stairs(total_events, ENERGY_GRID, color = "black", linewidth = 1.6, baseline = None, label = "Total LAr signal")
    plt.stairs(total_cryostat_events, ENERGY_GRID, label = "Cryostat contribution", baseline = None, color = "maroon", linewidth = 1.6)
    plt.axvline(5,color="gray",linestyle="--",linewidth=1.2,alpha=0.8)
    plt.xlim(0,50)
    plt.xlabel("Deposited energy in LAr [MeV]", fontsize = 12)
    plt.ylabel("Number of events per 0.5 MeV bin", fontsize = 12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.ylim(bottom=0)

    plt.legend(fontsize = 12)
    plt.text(0.47, 0.44, rf"""LAr events > 5 MeV: {total_events_mask.sum():.1f}
Smeared cryostat events: {total_cryostat_events.sum():.1f}
Cryostat events > 5 MeV: {total_cryostat_mask.sum():.1f}""",transform=plt.gca().transAxes, linespacing = 1.4, fontsize = 12)#,ha="right",va="top")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"{detector}, {flux}, {iron_thickness}mm, all iron channel")
    _savefig(f"plotting/{flux}_{detector}_argon_all_channels_and_cryostat_smeared.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

# Calculate expected event counts above 5 MeV for all four FD modules
def expected_event_counts(flux="garching", argon_detector="ar40kt"):

    mask = ENERGY_BIN_CENTERS > 5
    totals = {}

    for layer, info in DATA_LOOKUP.items():

        spectrum = np.zeros(len(ENERGY_BIN_CENTERS))

        for channel, channel_info in info["channels"].items():

            events = load_events_data(flux, get_detector(layer), channel, True)

            spectrum += events * channel_info["weight"]

        totals[layer] = 4 * spectrum[mask].sum()

    totals["Total cryostat"] = sum(totals.values())

    argon = sum(load_events_data(flux, argon_detector, channel, True) * weight for channel, weight in ARGON_CHANNELS.items())

    totals["Direct LAr reference"] = argon[mask].sum()

    return totals


def print_event_totals(flux="garching", argon_detector="ar40kt"):

    totals = expected_event_counts(flux, argon_detector)

    for name, count in totals.items():
        print(f"{name}: {count:.4f} events")


if __name__ == "__main__":

    print_event_totals("garching", "ar40kt")

    #plot_all_argon_channels_smeared_and_cryostat("garching", "ar40kt")
    #plot_all_iron_smeared_vs_unsmeared("garching", "PrimaryIron")
    #plot_one_channel_smeared_vs_unsmeared("garching", "nue_Fe56", "PrimaryIron")
    #plot_distribution_matrix("nue_Fe56", "PrimaryIron", False)
    #plot_deposition_efficiency("nue_Fe56", "PrimaryIron", 5)
    #plot_unsmeared_events_iron("garching", "FD_cryostat_primary")