"""Paths and unchanged scientific configuration for historical and final studies."""
import numpy as np
from pathlib import Path
import os

IRON_THICKNESS_MM = 62.7 #1.2mm lower bound
# Electron scattering is per atomic electron: Ar has 18, Fe has 26.
ARGON_CHANNELS = {
        "nue_Ar40" : 1,
        "nuebar_Ar40" : 1,
        "nue_e" : 18,
        "nuebar_e" : 18,
        "numu_e" : 18,
        "numubar_e" : 18,
        "nutau_e" : 18,
        "nutaubar_e" : 18,
        #"nc_nue_Ar40" : 1,
        #"nc_nuebar_Ar40" : 1,
        #"nc_numu_Ar40" : 1,
        #"nc_numubar_Ar40" : 1,
        #"nc_nutau_Ar40" : 1,
        #"nc_nutaubar_Ar40" : 1

    }



IRON_CHANNELS = {
        "nue_e" : 26,
        "nuebar_e" : 26,
        "numu_e" : 26,
        "numubar_e" : 26,
        "nutau_e" : 26,
        "nutaubar_e" : 26,
        "nue_Fe56" : 1,
        "nuebar_Fe56" : 1,
        "nc_nue_Fe56" : 1,
        "nc_nuebar_Fe56" : 1,
        "nc_numu_Fe56" : 1,
        "nc_numubar_Fe56" : 1,
        "nc_nutau_Fe56" : 1,
        "nc_nutaubar_Fe56" : 1,
    }

ENERGY_GRID = np.linspace(0.5, 100, 201)


PROJECT_ROOT = Path(__file__).resolve().parent

# Environment overrides preserve external simulation data in its existing location.
GEANT4_BUILD_DIR = Path(os.environ.get("GEANT4_BUILD_DIR", PROJECT_ROOT / "Geant4" / "build")).expanduser()
ROOT_DATA = Path(os.environ.get("FDSTEEL_ROOT_DATA", PROJECT_ROOT / "data" / "root")).expanduser()
SNOWGLOBES_DIR = Path(os.environ.get("SNOWGLOBES_DIR", PROJECT_ROOT / "snowglobes")).expanduser()
PLOT_DIR = Path(os.environ.get("PLOT_DIR", PROJECT_ROOT / "figures" / "python")).expanduser()
SNOWGLOBES = SNOWGLOBES_DIR  # Compatibility for existing Python callers.
ENERGY_BIN_CENTERS = 0.5 * (ENERGY_GRID[:-1] + ENERGY_GRID[1:])

# Historical only-iron model; thicknesses are in mm.
# Neutral-current entries in IRON_CHANNELS are historical, not final simulations.
# Muon/tau electron scattering shares the same neutral-current response at
# these energies; antineutrinos likewise reuse the corresponding antimuon files.
#simplified model only one iron layer 
MODEL_INFORMATION = {
    "nue_Fe56": {
        "pdg": 12,
        "scan_min_energy": 6.31,
        "weight": 1,

        "simulations": {
            12: {       #in mm
                "folder": ROOT_DATA /"12mm"/ "marley_CC_nue_Fe56_100k_12mm",
                "filename": "MARLEY_{E:.5f}_CC_nue_ALLOWED_Fe56_12mm.root",
            },

            1.2:
            {
                "folder": ROOT_DATA/"1.2mm"/"marley_Fe56_100k_1.2mm_nue_Fe56",
                "filename": "MARLEY_{E:.5f}_nue_Fe56_ALLOWED_Fe56_1.2mm.root",
            },
            62.7:
            {
                "folder": ROOT_DATA/"62.7mm"/"marley_Fe56_100k_62.7mm_nue_Fe56",
                "filename": "MARLEY_{E:.5f}_nue_Fe56_ALLOWED_Fe56_62.7mm.root",
            }
        },
    },

    "nuebar_Fe56": {
        "pdg": -12,
        "scan_min_energy": 5.1,
        "weight": 1,

        "simulations": {
            12: 
            {     #in mm
                "folder": ROOT_DATA / "12mm"/"marley_CC_nuebar_Fe56_100k_12mm",
                "filename": "MARLEY_{E:.5f}_CC_nuebar_ALLOWED_Fe56_12mm.root",
            },

            1.2:
            {
                "folder": ROOT_DATA/"1.2mm"/"marley_Fe56_100k_1.2mm_nuebar_Fe56",
                "filename": "MARLEY_{E:.5f}_nuebar_Fe56_ALLOWED_Fe56_1.2mm.root",
            },
            62.7:
            {
                "folder": ROOT_DATA/"62.7mm"/"marley_Fe56_100k_62.7mm_nuebar_Fe56",
                "filename": "MARLEY_{E:.5f}_nuebar_Fe56_ALLOWED_Fe56_62.7mm.root",
            }            
        },
    },
    "nue_e":{
        "pdg": 12,
        "scan_min_energy": 0.5,
        "weight": 26,

        "simulations":{
            12:{
                "folder": ROOT_DATA /"12mm"/ "marley_Fe56_100k_12mm_nue_e",
                "filename": "MARLEY_{E:.5f}_nue_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_nue_e",
                "filename": "MARLEY_{E:.5f}_nue_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:
            {
            "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_nue_e",
                "filename": "MARLEY_{E:.5f}_nue_e_ALLOWED_Fe56_62.7mm.root",      
            }
        },

    },
    "numu_e":{
        "pdg": 14,
        "scan_min_energy": 0.5,
        "weight":26,
        
        "simulations":{
            12:{
                "folder": ROOT_DATA /"12mm"/ "marley_Fe56_100k_12mm_numu_e",
                "filename": "MARLEY_{E:.5f}_numu_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_numu_e",
                "filename": "MARLEY_{E:.5f}_numu_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:
            {
            "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_numu_e",
            "filename": "MARLEY_{E:.5f}_numu_e_ALLOWED_Fe56_62.7mm.root",                   
            }                        
        },

    },
    "nutau_e":{
        "pdg": 16,
        "scan_min_energy": 0.5,
        "weight": 26,

        "simulations":{
            12:{
                "folder": ROOT_DATA / "12mm"/"marley_Fe56_100k_12mm_nutau_e",
                "filename": "MARLEY_{E:.5f}_nutau_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_numu_e", #numu <--> nutau
                "filename": "MARLEY_{E:.5f}_numu_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:
            {
                "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_numu_e", #numu <--> nutau
                "filename": "MARLEY_{E:.5f}_numu_e_ALLOWED_Fe56_62.7mm.root",                    
            }                        
        },

    },
    "nuebar_e":{
        "pdg": -12,
        "scan_min_energy": 0.5,
        "weight": 26,

        "simulations":{
            12:{
                "folder": ROOT_DATA /"12mm"/ "marley_Fe56_100k_12mm_nuebar_e",
                "filename": "MARLEY_{E:.5f}_nuebar_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_nuebar_e",
                "filename": "MARLEY_{E:.5f}_nuebar_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:{
                "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_nuebar_e",
                "filename": "MARLEY_{E:.5f}_nuebar_e_ALLOWED_Fe56_62.7mm.root",                
            }            
        },

    },
    "numubar_e":{
        "pdg":-14,
        "scan_min_energy": 0.5,
        "weight": 26,

        "simulations":{
            12:{
                "folder": ROOT_DATA / "12mm"/"marley_Fe56_100k_12mm_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:{
                "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_ALLOWED_Fe56_62.7mm.root",                
            }                              
        },

    },
    "nutaubar_e":{
        "pdg": -16,
        "scan_min_energy": 0.5,
        "weight":26,

        "simulations":{
            12:{
                "folder": ROOT_DATA /"12mm"/ "marley_Fe56_100k_12mm_nutaubar_e",
                "filename": "MARLEY_{E:.5f}_nutaubar_e_ALLOWED_Fe56_12mm.root",
            },
            1.2:{
                "folder": ROOT_DATA /"1.2mm"/ "marley_Fe56_100k_1.2mm_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_ALLOWED_Fe56_1.2mm.root",                
            },
            62.7:{
                "folder": ROOT_DATA /"62.7mm"/ "marley_Fe56_100k_62.7mm_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_ALLOWED_Fe56_62.7mm.root",                
            }                             
        },

    },
}
# Final cryostat: LAr -> PrimaryIron -> PU -> OuterIron (thicknesses in mm).
#numu <-> nutau
#numubar <-> nutaubar
DATA_LOOKUP = {
    "PrimaryIron":{
        "thickness": 1.2,

        "channels": 
        {
            "nue_Fe56":
            {
                "scan_min_energy": 6.31,
                "pdg": 12,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_nue_Fe56",
                "filename": "MARLEY_{E:.5f}_nue_Fe56_PrimaryIron.root", 
                "weight":1,
            },
            "nuebar_Fe56":
            {
                "scan_min_energy": 5.1,
                "pdg": -12,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_nuebar_Fe56",
                "filename": "MARLEY_{E:.5f}_nuebar_Fe56_PrimaryIron.root", 
                "weight":1,
            },
            "nue_e":
            {
                "scan_min_energy": 0.5,
                "pdg":12,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_nue_e",
                "filename": "MARLEY_{E:.5f}_nue_e_PrimaryIron.root",
                "weight":26,
            },
            "nuebar_e":
            {
                "scan_min_energy": 0.5,
                "pdg":-12,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_nuebar_e",
                "filename": "MARLEY_{E:.5f}_nuebar_e_PrimaryIron.root",
                "weight":26,
            },
            "numu_e":
            {
                "scan_min_energy": 0.5,
                "pdg":14,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_numu_e",
                "filename": "MARLEY_{E:.5f}_numu_e_PrimaryIron.root",
                "weight":26,
            },            
            "numubar_e":
            {
                "scan_min_energy": 0.5,
                "pdg":-14,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_PrimaryIron.root",
                "weight":26,
            },            
            "nutau_e":
            {
                "scan_min_energy": 0.5,
                "pdg":16,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_numu_e",
                "filename": "MARLEY_{E:.5f}_numu_e_PrimaryIron.root",
                "weight":26,
            },        
            "nutaubar_e":
            {
                "scan_min_energy": 0.5,
                "pdg":-16,
                "folder": ROOT_DATA/"iron_PU"/"PrimaryIron"/"marley_scan100k_PrimaryIron_numubar_e",
                "filename": "MARLEY_{E:.5f}_numubar_e_PrimaryIron.root",
                "weight":26,
            },                
        },
    },
    "PU":
    {
        "thickness": 800,
        "channels":
        {
            "ibd":
            {
                "pdg": -12,
                "scan_min_energy": 1.804,
                "folder":ROOT_DATA/"iron_PU"/"PU"/"marley_scan100k_PU_ibd",
                "filename": "MARLEY_{E:.5f}_ibd_PU.root",
                "weight": 1,
            }
        }
    },
    "OuterIron":
    {
        "thickness": 54.63,
        "channels":
               {
                    "nue_Fe56":
                    {
                        "scan_min_energy": 6.31,
                        "pdg": 12,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_nue_Fe56",
                        "filename": "MARLEY_{E:.5f}_nue_Fe56_OuterIron.root",
                        "weight":1, 
                    },
                    "nuebar_Fe56":
                    {
                        "scan_min_energy": 5.1,
                        "pdg": -12,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_nuebar_Fe56",
                        "filename": "MARLEY_{E:.5f}_nuebar_Fe56_OuterIron.root", 
                        "weight":1,
                    },
                    "nue_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":12,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_nue_e",
                        "filename": "MARLEY_{E:.5f}_nue_e_OuterIron.root",
                        "weight":26,
                    },
                    "nuebar_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":-12,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_nuebar_e",
                        "filename": "MARLEY_{E:.5f}_nuebar_e_OuterIron.root",
                        "weight":26,
                    },
                    "numu_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":14,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_numu_e",
                        "filename": "MARLEY_{E:.5f}_numu_e_OuterIron.root",
                        "weight":26,
                    },            
                    "numubar_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":-14,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_numubar_e",
                        "filename": "MARLEY_{E:.5f}_numubar_e_OuterIron.root",
                        "weight":26,
                    },            
                    "nutau_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":16,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_numu_e",
                        "filename": "MARLEY_{E:.5f}_numu_e_OuterIron.root",
                        "weight":26,
                    },        
                    "nutaubar_e":
                    {
                        "scan_min_energy": 0.5,
                        "pdg":-16,
                        "folder": ROOT_DATA/"iron_PU"/"OuterIron"/"marley_scan100k_OuterIron_numubar_e",
                        "filename": "MARLEY_{E:.5f}_numubar_e_OuterIron.root",
                        "weight":26,
                    },                
                },

    }
}

if __name__ == "__main__":
    pass
