# DUNE Far Detector Cryostat Supernova Neutrino Study

This repository contains the simulation and analysis code developed to investigate the contribution of supernova neutrino interactions in the DUNE Far Detector cryostat to the observable signal in liquid argon.

## Overview

The study considers neutrino interactions in the materials surrounding the liquid argon volume of the DUNE Far Detector. A simplified cryostat geometry is modeled, consisting of three layers:

- **Primary iron membrane:** 1.2 mm
- **Polyurethane insulation:** 800 mm
- **Outer iron structure:** 54.63 mm effective thickness

The study combines three simulation tools:

- **SNOwGLoBES:** Calculation of expected supernova neutrino interaction rates in the different cryostat materials and liquid argon.
- **MARLEY:** Generation of low-energy neutrino interaction final states, including charged-current interactions on iron-56 and neutrino-electron elastic scattering.
- **Geant4:** Simulation of particle propagation and energy deposition in the liquid argon volume, using a simplified cryostat geometry.

A separate Geant4-based generator is used to approximate inverse beta decay interactions in the polyurethane layer.

The resulting detector response matrices are incorporated into SNOwGLoBES to estimate the contribution of cryostat interactions to the observable supernova neutrino signal.

## Main Results

For a supernova at 10 kpc using the Garching flux model, the study predicts approximately 15.8 additional events originating from neutrino interactions in the cryostat materials, compared to approximately 1047 direct liquid argon events, assuming an energy-deposition threshold of 5 MeV.

The cryostat contribution corresponds to approximately 1.5% of the direct liquid argon signal.

These results represent first-order estimates based on a simplified cryostat geometry, interaction models, and detector response.

A corresponding Technical Note describing the methodology and results in more detail is in preparation.

## Repository Structure

- `Geant4/advanced_model/` – Final cryostat geometry and particle propagation simulations (C++)
- `Geant4/only_iron/` – Historical simplified cryostat model containing only an iron layer (C++)
- `marley_files/` – MARLEY interaction configurations and reaction files
- `snowglobes_files/` – SNOwGLoBES detector and channel configurations
- `snowglobes_helper/` – Python scripts for generating detector response matrices
- `macros/` – Python scripts for generating Geant4 simulation macros
- `plotting/` – Python analysis and visualization
- `figures/` – Selected simulation results and plots
- `config.py` and `helper.py` – Shared simulation parameters and analysis functions

The repository contains the project-specific simulation and analysis code. The original simulations were performed using separate local installations of Geant4, MARLEY, and SNOwGLoBES.

The large ROOT simulation files are not included.

## Running the Simulations

The Geant4 models require external installations of Geant4, MARLEY, and HepMC3, together with CMake and a C++17-compatible compiler.

The Python analysis uses NumPy, Matplotlib, and uproot.

The simulation workflow consists of generating neutrino interactions, propagating the resulting particles through the cryostat geometry, constructing detector response matrices, and calculating the expected event rates with SNOwGLoBES.

Local paths in the simulation configurations and macros must be adapted to the corresponding software installations. Reproducing the complete analysis also requires the original ROOT simulation data.

## MARLEY Reaction Data

The charged-current electron-neutrino reaction file for iron-56 used in this study was provided by Sam Hedges and is not included in this repository.

The electron-antineutrino reaction file (`vebar56FeCC_Scott2014.react`) was prepared as part of this project using digitized experimental Gamow–Teller strength data from Scott et al. (2014).

The neutrino-electron elastic scattering reaction file (`ES.react`) originates from MARLEY.

## Third-Party Software

This project contains code adapted from Geant4 examples and modified code based on MARLEY.

The original copyright and license notices in these source files have been retained.

This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4).
