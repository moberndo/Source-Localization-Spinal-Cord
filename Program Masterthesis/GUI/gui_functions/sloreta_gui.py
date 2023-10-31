# ---- IMPORT PACKAGES ----
from pathlib import Path
import matplotlib.pyplot as plt
from os.path import isfile, isdir
from os import mkdir
from numpy import load

# ----- IMPORT CUSTOM FUNCTIONS -----
from .custom_functions import ModelSelection, LoadData, ProcessData
from .custom_functions import DipolesLeadfield, ElectrodePositions
from .custom_functions import sLORETA, ShowSolution, GetWaveformPoints

# ---- IMPORT PARAMETERS ---- 
from .params import params_elec, mat_conduc

def sloreta():

    user_ID = 'XYZ'
    num_dipoles_per_section = 2
    fs = 512
    alpha = 0.05


    mesh, params_geom, loc_C7, maxh_ = ModelSelection(user_ID, main_path)

    _ , data_SCP, data_marker, num_participants = LoadData(main_path,
                                                           user_folder)
    data_SCP, data_SCP_participants = ProcessData(data_SCP, data_marker,
                                                  num_participants, fs,
                                                  main_path, user_ID)
    
    electrodes, N_E = ElectrodePositions(loc_C7, params_elec, params_geom)

    # Create dipoles that are necessary for the Leadfield-Calculation
    dipoles, z_vals = DipolesLeadfield(params_geom, num_dipoles_per_section)
    # Solve the inverser problem
    J = sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles, electrodes, 
                user_ID, main_path, alpha)
    
    # Calculate the Waveform-Points, at which the solution is visualized
    P1, N1, P2 = GetWaveformPoints(J, data_SCP, fs)
    J_WP = [J[P1,:], J[N1,:], J[P2,:]]
    
    ShowSolution(J_WP, params_geom, z_vals, num_dipoles_per_section,
             view = 'frontal', save = True, path = images_folder)


    