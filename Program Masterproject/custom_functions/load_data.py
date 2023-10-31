""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Load Data Function
    Author: Markus E. Oberndorfer
####################################################################### """

from scipy.io import loadmat   
from os import chdir
from glob import glob


def LoadData(num_participant, main_path):
    folder_participant = 'participant_' + num_participant
    folder_data = num_participant + '_data'
    path_participant_data = main_path / 'data' / folder_participant / folder_data
    
    chdir(path_participant_data)
    data_EEG = []
    data_SCP = []
    for file in glob("*.mat"):
        #print('Found: ', file)
        data_EEG.append(loadmat(file)['EEGdata'][2:18])
        data_SCP.append(loadmat(file)['EEGdata'][18:34])

    print('Data loaded.')
    return data_EEG, data_SCP


if __name__ == '__main__':
    from pathlib import Path
    num_participant = 'XY'
    main_path = Path('/Users/markusoberndorfer/Uni/Masterarbeit/Program')
    folder_participant = 'participant_' + num_participant
    folder_data = num_participant + '_data'
    path_participant_data = main_path / 'data' / folder_participant / folder_data
    data_EEG, data_SCP = LoadData(path_participant_data)