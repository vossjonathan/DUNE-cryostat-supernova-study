import numpy as np
from pathlib import Path


IRON_THICKNESS_MM = 62.7 #1.2mm lower bound
#IRON_MASS_T = calculate_iron_mass(IRON_THICKNESS_MM) #in t
#SNOWGLOBES_IRON_MASS = 2300 #in t DIFFERENT SNOWGLOBES MODELS FOR ALL THICKNESSES
#MASS_SCALE = IRON_MASS_T / SNOWGLOBES_IRON_MASS
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

ROOT_DATA = Path.home()/"software"/"geant4"/"geant4-work"/"FDSteel"/"build"/"root_data"


MODEL_INFORMATION = {
    "nue_Fe56": {
        "pdg": 12,
        "scan_min_energy": 6.31,
        "weight": 1,

        "simulations": {
            12: {       #in mm
                "folder": ROOT_DATA /"12mm"/ "marley_CC_nue_Fe56_100k_12mm",
                "filename": "MARLEY_{E:.5f}_CC_nue_ALLOWED_Fe56_12mm.root",
                #"folder": ROOT_DATA / "marley_CC_nue_Fe56_10k_12mm_no_8",
                #"filename": "MARLEY_{E:.5f}_CC_nue_ALLOWED_Fe56_12mm.root",
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

if __name__ == "__main__":
    pass