"""####################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Main
    Author: Markus E. Oberndorfer
####################################################################"""

# ---- IMPORT PACKAGES ----
from pathlib import Path
import matplotlib.pyplot as plt
from os.path import isfile, isdir
from os import mkdir
from numpy import load

# ----- IMPORT CUSTOM FUNCTIONS -----
from custom_functions import ModelSelection, LoadData, ProcessData
from custom_functions import DipolesLeadfield, ElectrodePositions
from custom_functions import ChannelPlotWP_1, ChannelPlotWP_2
from custom_functions import GetWaveformPoints, sLORETA, ShowSolution
from custom_functions import CalcHyperparams, VisualizeHyperParams

# ---- IMPORT PARAMETERS ---- 
from params import params_elec, mat_conduc

######### SET THESE VARIABLES #########
#-------------------------------------#
user_ID = 'XYZ'
num_dipoles_per_section = 2
fs = 512
alpha = 0.05
# These parameters can be changed if needed
HYPERPARAM_SEARCH = True
PLOT_HP_DATA = True
#------------------------------------#

################## MAIN ###################
print('\n--------------------------------------------')
print('Starting "Source Localization: Spinal Cord"')
print('-------------------------------------------- \n \n')

''' Setup '''
# Create folder structure
main_path = Path(__file__).parent
user_folder = 'user_' + user_ID
user_folder = main_path / 'data' / user_folder
images_folder = user_folder / 'images'
processed_data_folder = user_folder / 'processed_data'
if not isdir(user_folder):
    exit('Create a user folder. See ReadMe.md file.')
if not isdir(images_folder):
    mkdir(images_folder)
if not isdir(processed_data_folder):
    mkdir(processed_data_folder)

# Make paths for processed data
processed_data_path = main_path / 'data' / user_folder / 'processed_data' 
processed_data_path = processed_data_path / 'data_avg.npy'

processed_data_path_par = main_path / 'data' / user_folder / 'processed_data'
processed_data_path_par = processed_data_path_par / 'data_participants.npy'

#########################################
'''' ======== LOAD ALL FILES ======== '''
#########################################

# Select the model or load existing model
mesh, params_geom, loc_C7, maxh_ = ModelSelection(user_ID, main_path)

# Load the data and process it or load already processed data
if isfile(processed_data_path) and isfile(processed_data_path_par):
    data_SCP = load(processed_data_path)
    data_SCP_participants = load(processed_data_path_par)
else:
    _ , data_SCP, data_marker, num_participants = LoadData(main_path,
                                                           user_folder)
    data_SCP, data_SCP_participants = ProcessData(data_SCP, data_marker,
                                                  num_participants, fs,
                                                  main_path, user_ID)

# Plot the average signal with the waveform points being marked
ChannelPlotWP_1(data_SCP, fs, save = True, main_path = main_path,
                 user_ID = user_ID)
ChannelPlotWP_2(data_SCP, fs, save = True, main_path = main_path,
                 user_ID = user_ID)

# Define the position of the electrodes (depending on C7)
electrodes, N_E = ElectrodePositions(loc_C7, params_elec, params_geom)

############################################
'''' ==== CALCULATE HYPERPARAMETERS ==== '''
############################################

if HYPERPARAM_SEARCH:
    # These following two lists are the hyperparameter values that are tested
    alpha_hp = [0.005, 0.01, 0.05, 0.1, 0.2, 0.35, 0.5]
    num_dip_hp = [1,2,3,4,5]

    if PLOT_HP_DATA:
        hp, data_plot_hp = CalcHyperparams(params_geom, num_dip_hp, alpha_hp,
                                           electrodes, user_ID, main_path,
                                           data_SCP, mesh, mat_conduc, maxh_,
                                           fs, plot_data = True)
        F_vals_P1 = data_plot_hp[0]
        F_vals_N1 = data_plot_hp[1]
        F_vals_P2 = data_plot_hp[2]
        
        # Visualization of the Hyperparameters' Influence on
        # the Error-Functional
        VisualizeHyperParams(F_vals_P1, num_dip_hp, alpha_hp, 'P1',
                              save=True, main_path=main_path,
                              user_ID = user_ID)
        VisualizeHyperParams(F_vals_N1, num_dip_hp, alpha_hp, 'N1',
                              save=True, main_path=main_path,
                              user_ID = user_ID)
        VisualizeHyperParams(F_vals_P2, num_dip_hp, alpha_hp, 'P2',
                              save=True, main_path=main_path,
                              user_ID = user_ID)
    else:
        hp = CalcHyperparams(params_geom, num_dip_hp, alpha_hp,
                            electrodes, user_ID, main_path,
                            data_SCP, mesh, mat_conduc, maxh_,
                            fs, plot_data = False)
    
    alpha = hp[0]
    num_dipoles_per_section = hp[1]

#########################################
'''' ==== SOLVE INVERSE PROBLEM ==== '''
#########################################

# Create dipoles that are necessary for the Leadfield-Calculation
dipoles, z_vals = DipolesLeadfield(params_geom, num_dipoles_per_section)
# Solve the inverser problem
J = sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles, electrodes, 
            user_ID, main_path, alpha)

print('Inverse Problem calculated.')

#####################################
'''' ==== VISUALIZE SOLUTION ==== '''
#####################################

# Calculate the Waveform-Points, at which the solution is visualized
P1, N1, P2 = GetWaveformPoints(J, data_SCP, fs)
J_WP = [J[P1,:], J[N1,:], J[P2,:]]

# Create a 3D image of the solution at all time points P1, N1 and P2
# view options:
# frontal, lower_frontal, upper_frontal, upper_right, upper_left 
ShowSolution(J_WP, params_geom, z_vals, num_dipoles_per_section,
             view = 'frontal', save = True, path = images_folder)

plt.show()

print('\n--------------------------------------------')
print('Finished "Source Localization: Spinal Cord"')
print('-------------------------------------------- \n \n')
