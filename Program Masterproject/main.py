""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Main
    Author: Markus E. Oberndorfer
####################################################################### """

################# IMPORTS #################
''' ===== IMPORT Packages ===== '''
from sys import exit
from pathlib import Path
import matplotlib.pyplot as plt

''' ===== IMPORT Custom Functions ===== '''
from custom_functions import ModelSelection, LoadData, ProcessData, CreateArtificialSignal
from custom_functions import DipolePositions, ElectrodePositions
from custom_functions import sLORETA, CreateAnimation

''' ==== IMPORT params ===='''
from params import params_elec, mat_conduc

################## MAIN ###################
print('\n--------------------------------------------')
print('Starting "Source Localization: Spinal Cord"')
print('-------------------------------------------- \n \n')

'''
inv_method = input('Which method should be used to solve the inverse problem? \nType "Help" for more information. \n')
num_participant = input('Enter participant number (two digits) \n')
electrode_pattern = input('Which electrode pattern should be used? [16/32] \n')
num_dipoles = input('How many dipoles should be used? \n')
'''

inv_method = 'sLORETA'
num_participant = 'XY'
num_dipoles = 71
electrode_pattern = '16'

#########################################
'''' ======== LOAD ALL FILES ======== '''
#########################################

''' Select the model '''
main_path = Path(__file__).parent
mesh, params_geom, loc_C7, maxh_ = ModelSelection(num_participant, main_path, show = True)
print(loc_C7)

''' Load the data and process it '''
_ , data_SCP = LoadData(num_participant, main_path)  # first output is the EEG data
data_SCP = ProcessData(data_SCP)    # data_SPC is a list in which each entry contains one run
data_SCP = data_SCP[4]              # for now work with run 4 of

#########################################
'''' ==== CALCULATE DIPOLE AND ELECTRODE SETUP ==== '''
#########################################

''' Define how many dipoles should be used and where they are located '''
dipoles, num_dipoles, z_vals = DipolePositions(params_geom, num_dipoles)

''' Define the position of the electrodes (depending on C7) '''
electrodes, N_E = ElectrodePositions(electrode_pattern, loc_C7, params_elec, params_geom)

#########################################
'''' ==== CREATE ARTIFICIAL DATA ==== '''
#########################################

signal_length = 1000
dipole_loc = [[loc_C7]]
data_SCP_arti = CreateArtificialSignal(dipole_loc, electrodes, mesh, mat_conduc, maxh_,
                                       params_geom, signal_length, main_path, SNR = None)


#########################################
'''' ==== SOLVE INVERSE PROBLEM ==== '''
#########################################

if inv_method == 'sLORETA':
    J = sLORETA(data_SCP_arti, mesh, mat_conduc, maxh_, dipoles, electrodes, num_participant, main_path)
elif inv_method == 'DNN':   # Deep Neural Network
    exit('DNN not yet implemented. Program closes now.')
    # J = DNN(...)

#########################################
'''' ==== VISUALIZE SOLUTION ==== '''
#########################################

animator = CreateAnimation(J, params_geom, z_vals)

plt.show()

print('\n--------------------------------------------')
print('Finished "Source Localization: Spinal Cord"')
print('-------------------------------------------- \n \n')

'''
Classifier for left or right hand movement

'''