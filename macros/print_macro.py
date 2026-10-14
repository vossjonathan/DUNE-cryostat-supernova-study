from config import ENERGY_GRID
from pathlib import Path
import numpy as np
from config import MODEL_INFORMATION

FDSteel = Path.home()/"software"/"geant4"/"geant4-work"/"FDSteel"



def macro_Geant4_Marley(n_events, n_jobs, channel, iron_thickness):

    min_energy = MODEL_INFORMATION[channel]["scan_min_energy"]

    macro_dir = f"marley_scan_100k_{iron_thickness}mm_Fe_{channel}"
    output_dir = f"marley_Fe56_100k_{iron_thickness}mm_{channel}"

    bin_edges = ENERGY_GRID
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_centers = bin_centers[bin_centers >= min_energy]

    build_dir = FDSteel / "build"

    macro_path = build_dir / "macros"/macro_dir
    output_path = build_dir /"root_data"/f"{iron_thickness}mm"/output_dir



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
                    f"MARLEY_{E:.5f}_{channel}_ALLOWED_Fe56_{iron_thickness}mm.root",
                    file=f
                )

                print(f"/run/beamOn {n_events}\n", file=f)
        
    print("for i in {1..8}; do\n"
    f"./exampleB1 macros/{macro_dir}/marleyScan100k${{i}}.mac &\n"
    "done")



def macro_Geant4_Marley_missing(n_events=100000, n_jobs=8):

    bin_edges = np.asarray(ENERGY_GRID, dtype=float)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    build_dir = FDSteel / "build"

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

    # Redistribute missing energies over 8 NEW jobs
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
                    "root_data/marley_CC_nue_Fe56_100k_12mm/"
                    f"MARLEY_{E:.5f}_CC_nue_ALLOWED_Fe56_12mm.root",
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
    macro_Geant4_Marley(100000, 8, channel = "numubar_e", iron_thickness=62.7)
    #macro_Geant4_Marley_missing()