import matplotlib.pyplot as plt
import numpy as np
from config import IRON_CHANNELS, ENERGY_GRID, MODEL_INFORMATION, ARGON_CHANNELS
from helper import calculate_iron_mass
from pathlib import Path
import uproot
from matplotlib.colors import LogNorm

OUT = Path.home()/"software"/"snowglobes"/"out"
ROOT_DATA = Path.home()/"software"/"geant4"/"geant4-work"/"FDSteel"/"build"/"root_data"
SNOWGLOBES = Path.home()/"software"/"snowglobes"

#HELPER FUNCTION


def load_events_data(flux, detector, channel, smeared, suffix = "" ):
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

    return Edep

#ANALYSIS FUNCTIONS

def get_single_energy_distribution(channel, iron_thickness, energy, dep_edges):
    Edep = load_Edep(channel, iron_thickness, energy)
    n_generated = len(Edep)

    counts, _ = np.histogram(Edep, bins = dep_edges)
    column = counts  / n_generated
    return column

    
def build_distribution_matrix(channel, iron_thickness):
    channel_info = MODEL_INFORMATION[channel]
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
    plt.savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
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
    plt.savefig(f"plotting/{flux}_{detector}_unsmeared_spectrum.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def plot_distribution_matrix(channel, iron_thickness):

    matrix, dep_edges = build_distribution_matrix(channel, iron_thickness)
    print(matrix.sum(axis = 0))
    plt.imshow(matrix, origin = "lower", aspect = "auto",    extent=[ENERGY_GRID[0],ENERGY_GRID[-1],dep_edges[0],dep_edges[-1]], norm=LogNorm(vmin=1e-5, vmax=matrix.max()))
    plt.xlabel(r"True neutrino energy $E_\nu$ [MeV]")
    plt.ylabel(r"LAr deposited energy [MeV]")
    plt.title(f"LAr energy deposition distribution for {channel}, {iron_thickness}mm")
    plt.colorbar(label = "Probability")
    plt.savefig(f"plotting/energy_distribution_matrix_{channel}_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
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
    plt.savefig(f"plotting/{energy}_MeV_distribution_{channel}_{iron_thickness}mm.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def calculate_efficiencies(channel, iron_thickness, dep_threshold):
    scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]

    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

    efficiencies = np.zeros(len(energies))

    for i, E in enumerate(energies):
        if E < scan_min_energy:
            continue

        Edep = load_Edep(channel, iron_thickness, E)

        efficiencies[i] = np.mean(Edep > dep_threshold)

    return efficiencies


def plot_deposition_efficiency(channel, iron_thickness, dep_threshold):

    scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]
    energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    mask = energies >= scan_min_energy
    efficiencies = calculate_efficiencies(channel, iron_thickness, dep_threshold)

    plt.plot(energies[mask], efficiencies[mask])
    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.grid(True, alpha = 0.3)
    plt.title(f"Event fraction with deposition energy > {dep_threshold} MeV, {channel}, {iron_thickness}mm")
    plt.savefig(f"plotting/effiency_fraction_{channel}_{iron_thickness}mm_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
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
    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
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
    plt.savefig(f"plotting/{flux}_{channel}_{detector}_smeared_unsmeared_{iron_thickness}mm_5MeV.svg", bbox_inches="tight", dpi=900)
    plt.close()

def plot_all_iron_channels_smeared(flux, iron_thickness):
    detector = f"FDSteel_{iron_thickness}mm"

    total_events = np.zeros(200)
    for channel in MODEL_INFORMATION.keys():
        events_smeared = load_events_data(flux, detector, channel, smeared = True)
        weight = MODEL_INFORMATION[channel]["weight"]
        events_smeared_w = events_smeared *weight
        total_events += events_smeared_w
        plt.stairs(events_smeared_w, ENERGY_GRID, label = channel, baseline= None)

    print("unmasked: ", total_events)

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
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
    plt.savefig(f"plotting/{flux}_{detector}_iron_all_channels_smeared.svg", bbox_inches = "tight", dpi = 900)
    plt.close()
    total_argon_events = 1047
    print("There are a total of ", total_argon_events, " events with deposited energy > 5 MeV on argon.")
    print("For the cryostat model with an iron thickness of ", iron_thickness, "mm, there is a total of " \
    fr"{total_events.sum():.2f} events. {total_events_mask.sum():.2f} events with deposited energy > 5 MeV.")
    frac_thres = total_events_mask.sum() / total_argon_events
    print(r"For deposited energy > 5 MeV, this is a contribution of " + f"{(frac_thres*100*4):.2f}%.")


def plot_all_argon_channels_smeared(flux, detector):

    total_events = np.zeros(200)
    for channel in ARGON_CHANNELS.keys():
        events_smeared = load_events_data(flux, detector, channel, smeared = True)
        weight = ARGON_CHANNELS[channel]
        events_smeared_w = events_smeared *weight
        total_events += events_smeared_w
        plt.stairs(events_smeared_w, ENERGY_GRID, label = channel, baseline= None)

    print("unmasked: ", total_events)

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
    mask = bin_centers > 5
    total_events_mask = total_events[mask]
    print("masked: ", total_events_mask)

    plt.stairs(total_events, ENERGY_GRID, color = "black", linewidth = 1.6, baseline = None, label = "Total")
    plt.xlim(0,60)
    plt.xlabel("Deposited energy in LAr [MeV]")
    plt.ylabel("Number of events per 0.5 MeV bin")
    plt.legend()
    plt.text(0.52, 0.36, rf"""Total events smeared: {total_events.sum():.2f}
Events with $E_{{th}}$ = 5 MeV: {total_events_mask.sum():.2f}""",transform=plt.gca().transAxes, linespacing = 1.4)#,ha="right",va="top")
    plt.grid(True, alpha = 0.3)
    #plt.title(f"{detector}, {flux}, {iron_thickness}mm, all iron channel")
    plt.savefig(f"plotting/{flux}_{detector}_argon_all_channels_smeared.svg", bbox_inches = "tight", dpi = 900)
    plt.close()

def compare_deposition_effiencies(channel, dep_threshold):
    iron_thicknesses = [1.2, 12, 62.7]
    print("Threshold: ", dep_threshold)
    for iron_thickness in iron_thicknesses:
        print("Calculating for: ", iron_thickness)
        scan_min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]
        energies = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
        mask = energies >= scan_min_energy
        efficiencies = calculate_efficiencies(channel, iron_thickness, dep_threshold)

        plt.plot(energies[mask], efficiencies[mask], label = f"{iron_thickness}mm")

    plt.xlabel("True neutrino energy [MeV]")
    plt.ylabel(f"Fraction of events with Edep > {dep_threshold} MeV")
    plt.legend()
    plt.grid(True, alpha = 0.3)
    plt.title(f"Comparison of event fraction with deposition energy > {dep_threshold} MeV, {channel}")
    plt.savefig(f"plotting/comparison_effiency_fraction_{channel}_{dep_threshold}_MeV.svg", bbox_inches = "tight", dpi = 900)
    plt.close()
    

def bar_plot_event_count_comparison(flux, channel):

    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
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

    plt.savefig(f"plotting/{flux}_{channel}_event_count_comparison.svg",bbox_inches="tight",dpi=900)

    plt.close()

def bar_plot_all_channels(flux):
    bin_centers = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])
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

    plt.savefig(f"plotting/{flux}_all_channels_event_count_comparison.svg",bbox_inches="tight",dpi=900)

    plt.close()


if __name__ == "__main__":
    #unsmeared_events("garching", "FD_modul_Fe")
    #print(load_events_data("garching", "FD_modul_Fe", "nue_e", False))
    #energy = ENERGY_GRID[0]
    #plot_one_channel_smeared_unsmeared("garching", "nue_Fe56", "FDSteel_62.7mm", 62.7)
    '''CHANNEL = "nutaubar_e"
    plot_distribution_matrix(CHANNEL, 62.7)
    plot_deposition_efficiency(CHANNEL, iron_thickness=62.7, dep_threshold=0.0)
    plot_deposition_efficiency(CHANNEL, iron_thickness=62.7, dep_threshold=5.0)
    #plot_unsmeared_events_iron("garching", "FDSteel_1.2mm", 1.2)'''

    bar_plot_all_channels("garching")
    #bar_plot_event_count_comparison("garching", "nue_Fe56")
    #plot_all_iron_channels_smeared("garching", 62.7)
    #compare_deposition_effiencies("nue_Fe56", 5)
    #compare_deposition_effiencies("nue_Fe56", 0)
    #compare_deposition_effiencies("nue_Fe56", 1 )
    #compare_deposition_effiencies("nuebar_Fe56", 5)
    #compare_deposition_effiencies("nuebar_Fe56", 0)
    #compare_deposition_effiencies("nuebar_Fe56", 1 )
    #plot_unsmeared_events_argon("garching", "ar40kt")
    #plot_all_argon_channels_smeared("garching", "ar40kt")
    #plot_all_iron_channels_smeared("garching", 1.2)
    #for channel in MODEL_INFORMATION.keys():
     #   calculate_visible_events(channel, IRON_THICKNESS_MM, 5.0)
    #plot_unsmeared_events("garching", "FD_modul_Fe", IRON_THICKNESS_MM)
    #CHANNEL = "nuebar_Fe56"
    #plot_distribution_matrix(CHANNEL, IRON_THICKNESS_MM)
    #plot_hisneutrino_energy(CHANNEL, IRON_THICKNESS_MM, 20.15125)
    #plot_deposition_efficiency("nue_Fe56", iron_thickness= 1.2, dep_threshold= 0.0)
    #plot_deposition_efficiency(CHANNEL, IRON_THICKNESS_MM, 1.0)
    #plot_deposition_efficiency("nue_Fe56", iron_thickness=1.2, dep_threshold=5.0)
    #calculate_visible_events(CHANNEL1, IRON_THICKNESS_MM, 0.0)
    #calculate_visible_events(CHANNEL1, IRON_THICKNESS_MM, 5.0)
    #calculate_visible_events(CHANNEL2, IRON_THICKNESS_MM, 0.0)
    #calculate_visible_events(CHANNEL2, IRON_THICKNESS_MM, 5.0)
    #channels = ["nue_e", "nuebar_e", "numu_e", "numubar_e", "nutau_e", "nutaubar_e"]
    #for channel in channels:
        #plot_different_efficiencies(channel, IRON_THICKNESS_MM)
     #   plot_distribution_matrix(channel, IRON_THICKNESS_MM)
    #plot_distribution_matrix("nue_e", IRON_THICKNESS_MM)
    #plot_hist_neutrino_energy("nue_e", IRON_THICKNESS_MM, 99.75125)
    #plot_deposition_efficiency("nue_e", IRON_THICKNESS_MM, 0.0)