from pathlib import Path
import numpy as np
from config import MODEL_INFORMATION, DATA_LOOKUP, GEANT4_BUILD_DIR, ENERGY_BIN_CENTERS

def _validate(n_events, n_jobs):
    for name, value in (("n_events", n_events), ("n_jobs", n_jobs)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value <= 0:
            raise ValueError(f"{name} must be a positive integer")


def macro_Geant4_Marley_only_iron(n_events, n_jobs, channel, iron_thickness, output_dir=None):
    _validate(n_events, n_jobs)
    if channel not in MODEL_INFORMATION or iron_thickness not in MODEL_INFORMATION[channel]["simulations"]:
        raise ValueError("Unsupported historical channel or iron thickness")
    simulation = MODEL_INFORMATION[channel]["simulations"][iron_thickness]

    min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]
    macro_dir = f"marley_scan_100k_{iron_thickness}mm_Fe_{channel}"
    bin_centers = ENERGY_BIN_CENTERS
    bin_centers = bin_centers[bin_centers >= min_energy]
    build_dir = GEANT4_BUILD_DIR
    macro_path = Path(output_dir) if output_dir is not None else build_dir / "macros"/macro_dir
    output_path = simulation["folder"]

    macro_path.mkdir(parents=True, exist_ok=True)
    output_path.mkdir(parents=True, exist_ok=True)

    for job in range(n_jobs):

        energies = bin_centers[job::n_jobs]

        with open(macro_path / f"marleyScan100k{job+1}.mac", "w") as f:

            print("/run/numberOfThreads 1", file=f)
            print("/run/initialize", file=f)
            print("/tracking/storeTrajectory 0", file=f)
            print("/run/printProgress 10000\n", file=f)

            for E in energies:

                print(f"/marley/nuEnergy {E:.5f} MeV", file=f)

                print(
                    f"/analysis/setFileName "
                    f"{output_path}/"
                    f"{simulation['filename'].format(E=E)}",
                    file=f
                )

                print(f"/run/beamOn {n_events}\n", file=f)

    print(f"for i in {{1..{n_jobs}}}; do\n"
    f'./exampleB1 "{macro_path}/marleyScan100k${{i}}.mac" &\n'
    "done")

def macro_Geant4_Marley_PU_iron(n_events, n_jobs, channel, layer, output_dir=None):
    _validate(n_events, n_jobs)
    if layer not in DATA_LOOKUP or channel not in DATA_LOOKUP[layer]["channels"]:
        raise ValueError("Unsupported final-model layer/channel")

    config = DATA_LOOKUP[layer]["channels"][channel]

    min_energy = config["scan_min_energy"]
    nuFlavor= config["pdg"]
    macro_dir = f"marley_scan100k_{layer}_{channel}"
    bin_centers = ENERGY_BIN_CENTERS
    bin_centers = bin_centers[bin_centers >= min_energy]
    build_dir = GEANT4_BUILD_DIR
    macro_path = build_dir / "macros"/macro_dir
    output_path = config["folder"]
    if output_dir is not None:
        macro_path = Path(output_dir)


    macro_path.mkdir(parents=True, exist_ok=True)
    output_path.mkdir(parents=True, exist_ok=True)

    for job in range(n_jobs):

        energies = bin_centers[job::n_jobs]

        with open(macro_path / f"marleyScan100k{job+1}.mac", "w") as f:

            print("/run/numberOfThreads 1", file=f)
            print("/run/initialize", file=f)
            print("/tracking/storeTrajectory 0", file=f)
            print("/run/printProgress 10000\n", file=f)
            print(f"/marley/ironVolume {layer}", file = f)
            print(f"/marley/nuFlavor {nuFlavor}", file = f)

            for E in energies:

                print(f"/marley/nuEnergy {E:.5f} MeV", file=f)

                print(
                    f"/analysis/setFileName "
                    f"{output_path}/"
                    f"{config['filename'].format(E=E)}",
                    file=f
                )

                print(f"/run/beamOn {n_events}\n", file=f)

    print(f"for i in {{1..{n_jobs}}}; do\n"
    f'./exampleB1 "{macro_path}/marleyScan100k${{i}}.mac" &\n'
    "done")

def macro_IBD(n_events, output_dir=None, n_jobs=8):
    _validate(n_events, n_jobs)

    config = DATA_LOOKUP["PU"]["channels"]["ibd"]
    min_energy = config["scan_min_energy"]

    bin_centers = ENERGY_BIN_CENTERS
    bin_centers = bin_centers[bin_centers >= min_energy]
    build_dir = GEANT4_BUILD_DIR
    macro_path = Path(output_dir) if output_dir is not None else build_dir / "macros"
    output_path = config["folder"]


    macro_path.mkdir(parents=True, exist_ok=True)
    output_path.mkdir(parents=True, exist_ok=True)

    with open(macro_path / f"IBD_scan.mac", "w") as f:

        print(f"/run/numberOfThreads {n_jobs}", file=f)
        print("/run/initialize", file=f)
        print("/tracking/storeTrajectory 0", file=f)
        print("/run/printProgress 10000\n", file=f)
        print(f"/FDSteel/channel IBD", file = f)

        for E in bin_centers:

            filename = config["filename"].format(E=E)

            print(f"/FDSteel/nuEnergy {E:.5f} MeV", file=f)

            print(
                f"/analysis/setFileName "
                f"{output_path}/"
                f"{filename}",
                file=f
            )

            print(f"/run/beamOn {n_events}\n", file=f)


def macro_Geant4_Marley_missing(n_events=100000, n_jobs=8):
    _validate(n_events, n_jobs)
    simulation = MODEL_INFORMATION["nue_Fe56"]["simulations"][12]
    simulation["folder"].mkdir(parents=True, exist_ok=True)

    bin_centers = ENERGY_BIN_CENTERS
    build_dir = GEANT4_BUILD_DIR
    macro_dir = (
        build_dir
        / "macros"
        / "marley_scan_100k_12mm_Fe_CC_nue_missing"
    )

    macro_dir.mkdir(parents=True, exist_ok=True)

    # Original scan used 8 jobs.
    original_n_jobs = 8

    # Original jobs 1-4 did not finish.
    missing_energies = []

    for original_job in range(4):

        energies = bin_centers[original_job::original_n_jobs]

        # Remove the low-energy points where MARLEY hangs
        energies = energies[energies >= 6.31]

        missing_energies.extend(energies)

    # Sort all missing energies
    missing_energies = np.sort(np.asarray(missing_energies))

    print(f"Total missing energies: {len(missing_energies)}")
    print(
        f"Range: {missing_energies[0]:.5f} "
        f"-> {missing_energies[-1]:.5f} MeV"
    )

    # Redistribute missing energies over the requested number of recovery jobs
    for job in range(n_jobs):

        energies = missing_energies[job::n_jobs]

        macro_file = macro_dir / f"marleyMissing{job+1}.mac"

        with open(macro_file, "w") as f:

            print("/run/numberOfThreads 1", file=f)
            print("/run/initialize", file=f)

            print("/tracking/storeTrajectory 0", file=f)
            print("/tracking/verbose 0", file=f)
            print("/run/verbose 0", file=f)
            print("/run/printProgress 10000", file=f)
            print("", file=f)

            for E in energies:

                print(
                    f"/marley/nuEnergy {E:.5f} MeV",
                    file=f
                )

                print(
                    "/analysis/setFileName "
                    f"{simulation['folder']}/{simulation['filename'].format(E=E)}",
                    file=f
                )

                print(
                    f"/run/beamOn {n_events}",
                    file=f
                )

                print("", file=f)

        if len(energies) > 0:
            print(
                f"Job {job+1}: "
                f"{len(energies)} energies, "
                f"{energies[0]:.5f} -> {energies[-1]:.5f} MeV"
            )


if __name__ == "__main__":

    # Generate IBD scan
    macro_IBD(100000)

    # Generate MARLEY scans for the final cryostat model
    #macro_Geant4_Marley_PU_iron(100000, 8, "nue_Fe56", "PrimaryIron")
    #macro_Geant4_Marley_PU_iron(100000, 8, "nue_Fe56", "OuterIron")

    # Historical only-iron model
    #macro_Geant4_Marley_only_iron(100000, 8, "nue_Fe56", 12)

    # Recover missing historical simulation points
    #macro_Geant4_Marley_missing()