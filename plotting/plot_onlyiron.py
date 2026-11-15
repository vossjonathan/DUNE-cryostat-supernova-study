import matplotlib.pyplot as plt
import numpy as np
from config import IRON_CHANNELS, ENERGY_GRID, MODEL_INFORMATION, ARGON_CHANNELS
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


#HELPER FUNCTIONS


def load_events_data(flux, detector, channel, smeared, suffix = "" ):

    ending = f"smeared_unweighted{suffix}" if smeared else "unweighted"
    file = OUT / f"{flux}_{channel}_{detector}_events_{ending}.dat"
    data = np.loadtxt(file, max_rows=len(ENERGY_BIN_CENTERS), ndmin=2)
    if data.shape[0] != len(ENERGY_BIN_CENTERS) or data.shape[1] < 2:
        raise ValueError(f"{file}: expected {len(ENERGY_BIN_CENTERS)} rows and at least two columns")
    if not np.isfinite(data).all() or (data[:, 1] < 0).any():
        raise ValueError(f"{file}: expected finite nonnegative event counts")
    return data[:, 1]


def get_root_file(channel, iron_thickness, energy):

    channel_info = MODEL_INFORMATION[channel]
    simulation_info = channel_info["simulations"][iron_thickness]

    folder = simulation_info["folder"]
    filename = simulation_info["filename"].format(E = energy)
    return folder / filename


def load_Edep(channel, iron_thickness, energy):

    root_file = get_root_file(channel, iron_thickness, energy)

    with uproot.open(root_file) as file:
        Edep = file["FDSteel"]["EdepLAr"].array(library="np")

    if Edep.ndim != 1 or not len(Edep) or not np.isfinite(Edep).all():
        raise ValueError(f"{root_file}: expected nonempty finite scalar deposited energies in MeV")
    return Edep

#ANALYSIS FUNCTIONS

def get_single_energy_distribution(channel, iron_thickness, energy, dep_edges):

    Edep = load_Edep(channel, iron_thickness, energy)
    n_generated = len(Edep)
    if not n_generated or not np.isfinite(Edep).all():
        raise ValueError("Expected nonempty, finite deposited energies in MeV")

    counts, _ = np.histogram(Edep, bins = dep_edges)
    column = counts  / n_generated
    return column


def build_distribution_matrix(channel, iron_thickness):

    channel_info = MODEL_INFORMATION[channel]
    scan_min_energy = channel_info["scan_min_energy"]

    energies = ENERGY_BIN_CENTERS
    dep_edges = ENERGY_GRID

    n_true_bins = len(energies)
    n_dep_bins = len(dep_edges) - 1
    distribution_matrix = np.zeros((n_dep_bins, n_true_bins))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue
        column = get_single_energy_distribution(channel, iron_thickness, E, dep_edges)
        distribution_matrix[:, i] = column

    return distribution_matrix, dep_edges

#PLOT FUNCTIONS


def plot_unsmeared_events_iron(flux, detector, iron_thickness):

    total_events = 0
    for channel in MODEL_INFORMATION.keys():
        events_load = load_events_data(flux, detector, channel, smeared = False)
        weight = MODEL_INFORMATION[channel]["weight"]
        events_weighted = events_load * weight
        total_events += sum(events_weighted)
        plt.stairs(events_weighted, ENERGY_GRID, label = f"{channel}")

    plt.legend()
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel("Events per 0.5 MeV bin")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"Unsmeared spectrum {detector}, {flux}, {iron_thickness}mm")
    plt.text(x = 0.71, y = 0.5, s = f"Total: {total_events:.2f}", transform = plt.gca().transAxes, ha = "left", va = "top")
    _savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_unsmeared_events_argon(flux, detector):

    total_events = 0
    for channel in ARGON_CHANNELS.keys():
        events_load = load_events_data(flux, detector, channel, smeared = False)
        weight = ARGON_CHANNELS[channel]
        events_weighted = events_load * weight
        total_events += sum(events_weighted)
        plt.stairs(events_weighted, ENERGY_GRID, label = f"{channel}")

    plt.legend()
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel("Events per 0.5 MeV bin")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"Unsmeared spectrum {detector}, {flux}, {iron_thickness}mm")
    plt.text(x = 0.71, y = 0.5, s = f"Total: {total_events:.2f}", transform = plt.gca().transAxes, ha = "left", va = "top")
    _savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_distribution_matrix(channel, iron_thickness):

    matrix, dep_edges = build_distribution_matrix(channel, iron_thickness)
    print(matrix.sum(axis = 0))
    plt.imshow(matrix, origin = "lower", aspect = "auto",    extent=[ENERGY_GRID[0],ENERGY_GRID[-1],dep_edges[0],dep_edges[-1]], norm=LogNorm(vmin=1e-5, vmax=matrix.max()))
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel(r"LAr deposited energy [MeV]")
    plt.title(f"LAr energy deposition distribution for {channel}, {iron_thickness}mm")
    plt.colorbar(label = "Probability")
    _savefig(f"plotting/energy_distribution_matrix_{channel}_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


def plot_hist_neutrino_energy(channel, iron_thickness, energy):

    Edep = load_Edep(channel, iron_thickness, energy)

    Ar_dep_energy_grid = np.arange(0, 100.5, 0.5)

    zero_fraction = np.mean(Edep == 0)
    positive_fraction = np.mean(Edep > 0)

    print(f"E_nu = {energy}, exact zero:", np.sum(Edep == 0))
    print(f"E_nu = {energy}, 0 < Edep < 0.5:", np.sum((Edep > 0) & (Edep < 0.5)))
    print(f"E_nu = {energy} Edep >= 0.5:", np.sum(Edep >= 0.5))

    print(f"Exactly zero: {zero_fraction:.4f}")
    print(f"Edep > 0:     {positive_fraction:.4f}")
    plt.hist(Edep, bins = Ar_dep_energy_grid)
    plt.title(fr"$E_{{dep}}$ in LAr distribution for $E_\nu$ = {energy} MeV, {channel}, {iron_thickness}mm")
    plt.xlabel("Edep [MeV]")
    plt.ylabel("Event counts")
    _savefig(f"plotting/{energy}_MeV_distribution_{channel}_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def calculate_efficiencies(channel, iron_thickness, dep_threshold):

    scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]

    energies = ENERGY_BIN_CENTERS

    efficiencies = np.zeros(len(energies))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue

        Edep = load_Edep(channel, iron_thickness, E)

        efficiencies[i] = np.mean(Edep > dep_threshold)

    return efficiencies


def plot_deposition_efficiency(channel, iron_thickness, dep_threshold):

    scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]
    energies = ENERGY_BIN_CENTERS
    mask = energies >= scan_min_energy
    efficiencies = calculate_efficiencies(channel, iron_thickness, dep_threshold)

    plt.plot(energies[mask], efficiencies[mask])
    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.grid(True, alpha = 0.3)
    plt.title(f"Event fraction with deposition energy > {dep_threshold} MeV, {channel}, {iron_thickness}mm")
    _savefig(f"plotting/effiency_fraction_{channel}_{iron_thickness}mm_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def calculate_visible_events(channel, iron_thickness, dep_threshold):

    detector = f"FDSteel_{iron_thickness}mm"
    events_load = load_events_data("garching", detector, channel, smeared = False)
    events_weighted = events_load * IRON_CHANNELS[channel]
    efficiencies = calculate_efficiencies(channel, iron_thickness, dep_threshold)
    events_visible = events_weighted * efficiencies
    print("Iron Thickness: ", iron_thickness, "mm")
    print("Edep threshold: ", dep_threshold)
    print("Channel: ", channel)
    print("Unsmeared events: ", sum(events_weighted))
    print("Visible events: ", sum(events_visible))
    print("Fraction: ", sum(events_visible) / sum(events_weighted))

    return sum(events_visible)

def plot_different_efficiencies(channel, iron_thickness):

    thresholds = [0.0, 0.5, 1.0, 5.0]
    print("Plotting for: ", channel)
    for threshold in thresholds:
        plot_deposition_efficiency(channel, iron_thickness, threshold)

def plot_one_channel_smeared_unsmeared(flux, channel, detector, iron_thickness):

    energy = ENERGY_GRID
    bin_centers = ENERGY_BIN_CENTERS
    weight = MODEL_INFORMATION[channel]["weight"]

    events_smeared = load_events_data(flux, detector, channel, smeared = True, suffix = "")
    events_unsmeared = load_events_data(flux, detector, channel, smeared = False)
    mask = bin_centers > 5

    events_smeared_w = events_smeared * weight
    events_smeared_5MeV_w = events_smeared[mask] * weight
    events_unsmeared_w = events_unsmeared * weight


    plt.stairs(events_smeared_w, energy, label="Smeared: LAr deposited E", color="blue")
    plt.stairs(events_unsmeared_w, energy, label=r"Unsmeared: true $E_\nu$", color="red")
    plt.xlim(0,50)
    plt.text(0.52, 0.7, rf"""Total events unsmeared: {events_unsmeared_w.sum():.2f}
Total events smeared: {events_smeared_w.sum():.2f}
Events with $E_{{dep}}$ > 5 MeV: {events_smeared_5MeV_w.sum():.2f}""",transform=plt.gca().transAxes, linespacing = 1.4)#,ha="right",va="top")
    #plt.text(0.84, 0.77,f"Total events unsmeared: {events_unsmeared.sum():.2f}",transform=plt.gca().transAxes,ha="right",va="top")
    plt.xlabel(r" E [MeV]")
    plt.ylabel("Number of Events per bin")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.title(f"Smeared vs. unsmeared signal, {flux} flux, {channel}, {detector}")
    _savefig(f"plotting/{flux}_{channel}_{detector}_smeared_unsmeared_{iron_thickness}mm_5MeV.svg", bbox_inches="tight", dpi=900)
    plt.close()

def plot_all_iron_channels_smeared(flux, iron_thickness):
    """Plot diagnostic or expected counts; energy arguments are in MeV."""
    detector = f"FDSteel_{iron_thickness}mm"

    total_events = np.zeros(len(ENERGY_BIN_CENTERS))
    for channel in MODEL_INFORMATION.keys():
        events_smeared = load_events_data(flux, detector, channel, smeared = True)
        weight = MODEL_INFORMATION[channel]["weight"]
        events_smeared_w = events_smeared *weight
        total_events += events_smeared_w
        plt.stairs(events_smeared_w, ENERGY_GRID, label = channel, baseline= None)

    print("unmasked: ", total_events)
    bin_centers = ENERGY_BIN_CENTERS
    mask = bin_centers > 5
    total_events_mask = total_events[mask]
    print("masked: ", total_events_mask)

    plt.stairs(total_events, ENERGY_GRID, color = "black", linewidth = 1.6, baseline = None, label = "Total")
    plt.xlim(0,20)
    plt.xlabel("Deposited energy in LAr [MeV]")
    plt.ylabel("Number of events per 0.5 MeV bin")
    plt.legend()
    plt.text(0.57, 0.36, rf"""Total events smeared: {total_events.sum():.2f}
Events with $E_{{th}}$ = 5 MeV: {total_events_mask.sum():.2f}""",transform=plt.gca().transAxes, linespacing = 1.4)#,ha="right",va="top")
    plt.grid(True, alpha = 0.3)
    plt.title(f"{detector}, {flux}, {iron_thickness}mm, all iron channel")
    _savefig(f"plotting/{flux}_{detector}_iron_all_channels_smeared.svg", bbox_inches = "tight", dpi = 900)
    plt.close()
    total_argon_events = calculate_argon_total(flux)
    print("There are a total of ", total_argon_events, " events with deposited energy > 5 MeV on argon.")
    print("For the cryostat model with an iron thickness of ", iron_thickness, "mm, there is a total of " \
    fr"{total_events.sum():.2f} events. {total_events_mask.sum():.2f} events with deposited energy > 5 MeV.")
    frac_thres = total_events_mask.sum() / total_argon_events
    print(r"For deposited energy > 5 MeV, this is a contribution of " + f"{(frac_thres*100*4):.2f}%.")


def plot_all_argon_channels_smeared(flux, detector):

    total_events = np.zeros(len(ENERGY_BIN_CENTERS))
    total_eES_events = np.zeros(len(ENERGY_BIN_CENTERS))

    for channel in ARGON_CHANNELS.keys():
        events_smeared = load_events_data(flux, detector, channel, smeared = True)
        weight = ARGON_CHANNELS[channel]
        events_smeared_w = events_smeared *weight
        if channel.endswith("e"):
            total_eES_events += events_smeared_w
        if channel == "nue_Ar40":
            plt.stairs(events_smeared_w, ENERGY_GRID, label = r"$\nu_e$ CC on Ar40", baseline= None)
        if channel == "nuebar_Ar40":
            plt.stairs(events_smeared_w, ENERGY_GRID, label = r"$\bar\nu_e$ CC on Ar40", baseline= None)
        total_events += events_smeared_w


    print("unmasked: ", total_events)
    bin_centers = ENERGY_BIN_CENTERS
    mask = bin_centers > 5
    total_events_mask = total_events[mask]
    print("masked: ", total_events_mask)


    plt.stairs(total_eES_events, ENERGY_GRID, label = "Total eES", baseline= None)

    plt.stairs(total_events, ENERGY_GRID, color = "black", linewidth = 1.6, baseline = None, label = "Total")
    plt.axvline(
    5,
    color="gray",
    linestyle="--",
    linewidth=1.2,
    alpha=0.8
)
    plt.xlim(0,50)
    plt.xlabel("Deposited energy in LAr [MeV]", fontsize = 12)
    plt.ylabel("Number of events per 0.5 MeV bin", fontsize = 12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.ylim(bottom=0)

    plt.legend(fontsize = 12)
    plt.text(0.48, 0.57, rf"""Total events smeared: {total_events.sum():.0f}
Events with $E_{{th}}$ = 5 MeV: {total_events_mask.sum():.0f}""",transform=plt.gca().transAxes, linespacing = 1.4, fontsize = 12)#,ha="right",va="top")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"{detector}, {flux}, {iron_thickness}mm, all iron channel")
    _savefig(f"plotting/{flux}_{detector}_argon_all_channels_smeared.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def compare_deposition_effiencies(channel, dep_threshold):

    iron_thicknesses = [1.2, 12, 62.7]
    print("Threshold: ", dep_threshold)
    for iron_thickness in iron_thicknesses:
        print("Calculating for: ", iron_thickness)
        scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]
        energies = ENERGY_BIN_CENTERS
        mask = energies >= scan_min_energy
        efficiencies = calculate_efficiencies(channel, iron_thickness, dep_threshold)

        plt.plot(energies[mask], efficiencies[mask], label = f"{iron_thickness}mm")

    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.legend()
    plt.grid(True, alpha = 0.3)
    plt.title(f"Comparison of event fraction with deposition energy > {dep_threshold} MeV, {channel}")
    _savefig(f"plotting/comparison_effiency_fraction_{channel}_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
    plt.close()


def bar_plot_event_count_comparison(flux, channel):

    bin_centers = ENERGY_BIN_CENTERS
    mask_5MeV = bin_centers > 5

    weight = MODEL_INFORMATION[channel]["weight"]

    thicknesses = [1.2, 12, 62.7]
    detectors = ["FDSteel_1.2mm","FDSteel_12mm","FDSteel_62.7mm"]

    raw_counts = []
    smeared_counts = []
    threshold_counts = []

    for detector in detectors:

        raw = load_events_data(flux,detector,channel,smeared=False)* weight


        smeared = load_events_data(flux,detector,channel,smeared=True)* weight

        raw_counts.append(raw.sum())
        smeared_counts.append(smeared.sum())
        threshold_counts.append(smeared[mask_5MeV].sum())

    x = np.arange(len(thicknesses))
    width = 0.25


    bars_raw = plt.bar(x - width,raw_counts,width,label="Unsmeared")

    bars_smeared = plt.bar(x,smeared_counts,width,label="Smeared")

    bars_threshold = plt.bar(x + width, threshold_counts, width, label=r"$E_{\rm dep}>5$ MeV")

    plt.xticks(x,[f"{t:g} mm" for t in thicknesses])

    plt.xlabel("Iron thickness")
    plt.ylabel("Expected number of events")

    plt.title(f"Event counts vs. iron thickness {channel}, {flux} flux")

    plt.legend()


    # Numbers above bars
    for bars in [bars_raw,bars_smeared,bars_threshold]:
        plt.bar_label(bars,fmt="%.2f",padding=3,fontsize=8)

    plt.tight_layout()

    _savefig(f"plotting/{flux}_{channel}_event_count_comparison.svg",bbox_inches="tight",dpi=900)

    plt.close()

def bar_plot_all_channels(flux):

    bin_centers = ENERGY_BIN_CENTERS
    mask_5MeV = bin_centers > 5


    thicknesses = [1.2, 12, 62.7]
    detectors = ["FDSteel_1.2mm","FDSteel_12mm","FDSteel_62.7mm"]

    raw_counts = []
    smeared_counts = []
    threshold_counts = []

    for detector in detectors:
        total_raw = 0
        total_smeared = 0
        total_thres = 0
        for channel in MODEL_INFORMATION.keys():
            weight = MODEL_INFORMATION[channel]["weight"]

            raw = load_events_data(flux,detector,channel,smeared=False)* weight
            smeared = load_events_data(flux,detector,channel,smeared=True)* weight
            thres = smeared[mask_5MeV]

            total_raw += raw.sum()
            total_smeared += smeared.sum()
            total_thres += thres.sum()

        raw_counts.append(total_raw)
        smeared_counts.append(total_smeared)
        threshold_counts.append(total_thres)

    x = np.arange(len(thicknesses))
    width = 0.25


    bars_raw = plt.bar(x - width,raw_counts,width,label="Raw interactions")

    bars_smeared = plt.bar(x,smeared_counts,width,label="Smeared")

    bars_threshold = plt.bar(x + width, threshold_counts, width, label=r"$E_{\rm dep}>5$ MeV")

    plt.xticks(x,[f"{t:g} mm" for t in thicknesses])
    plt.ylim(0, 1.2 * max(raw_counts))

    plt.xlabel("Iron thickness")
    plt.ylabel("Expected number of events")

    #plt.title(f"Event counts vs. iron thickness all channels, {flux} flux")

    plt.legend()
    for bars in [bars_raw,bars_smeared,bars_threshold]:
        plt.bar_label(bars,fmt="%.2f",padding=3,fontsize=8)

    plt.tight_layout()

    _savefig(f"plotting/{flux}_all_channels_event_count_comparison.svg",bbox_inches="tight",dpi=900)

    plt.close()


def calculate_argon_total(flux="garching", detector="ar40kt"):

    spectrum = sum(load_events_data(flux, detector, channel, True) * weight for channel, weight in ARGON_CHANNELS.items())

    return spectrum[ENERGY_BIN_CENTERS > 5].sum()


if __name__ == "__main__":

    plot_all_argon_channels_smeared("garching", "ar40kt")

    #plot_all_iron_channels_smeared("garching", 12)
    #plot_distribution_matrix("nue_Fe56", 12)
    #plot_one_channel_smeared_unsmeared("garching", "nue_Fe56", "FDSteel_12mm", 12)
    #plot_deposition_efficiency("nue_Fe56", 12, 5)
    #bar_plot_all_channels("garching")