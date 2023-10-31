""" #######################################################################
    Title: Construction of an FE model of the spinal cord to
    simulate the voltage distribution caused by a dipole
    Subtitle: Load Data Function
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from scipy.io import loadmat   
from os import chdir
from os.path import join
from glob import glob

# ----- FUNCTIONS  -----
def LoadData(main_path, user_folder): 
    """
    Loads the data from a given path.
    
    Input:
        - main_path: Contains the path of the main file [pathlib object]
    
    Output:
        - data_participants_EEG: EEG data of the experiment. [Numpy array]
        - data_participants_SCP: SPC data of the experiment. [Numpy array]
        - data_participants_marker: Marker data that was sent during the
        experiment. [Numpy array]
        - num_participants:  Number of participants that took part in the 
        experiment. [int]
    """  
    def LoadParticipantData(participant_folder_path):
        """
        Loads the data of a singel participant.

        Input:
            - participant_folder_path: Path to the participant folder that
            stores the data regarding this participant. [Pathlib object]
        
        Output: 
            - data_marker: Marker data that was sent during the
            experiment. [Numpy array]
            - data_EEG: EEG data of the experiment. [Numpy array]
            - data_SCP: SPC data of the experiment. [Numpy array]
        """
        chdir(participant_folder_path)
        data_marker = []
        data_EEG = []
        data_SCP = []
        for file in glob("*.mat"):
            data_marker.append(loadmat(file)['EEGdata'][0])
            data_EEG.append(loadmat(file)['EEGdata'][2:18])
            data_SCP.append(loadmat(file)['EEGdata'][18:34])
        return data_marker, data_EEG, data_SCP


    # save all data to these lists
    data_participants_marker = []
    data_participants_EEG = []
    data_participants_SCP = []
    # find all participant folders
    participant_folder_path = main_path / 'data' / user_folder / 'SpinalCordPotentials'
    participant_folders_path = glob(join(participant_folder_path, "*", ""))
    num_participants = len(participant_folders_path)
    print(f'{num_participants} participants found.')

    for par_folder_path in participant_folders_path:
        data_marker, data_EEG, data_SCP = LoadParticipantData(par_folder_path)
        data_participants_marker.append(data_marker)
        data_participants_EEG.append(data_EEG)
        data_participants_SCP.append(data_SCP)

    print('All data loaded.')
    return data_participants_EEG, data_participants_SCP, \
           data_participants_marker, num_participants


if __name__ == '__main__':
    from pathlib import Path
    num_participant = 'XY'
    main_path = Path('/Users/markusoberndorfer/Uni/Masterarbeit/Program')
    folder_participant = 'participant_' + num_participant
    folder_data = num_participant + '_data'
    path_participant_data = main_path / 'data' / folder_participant / folder_data
    data_EEG, data_SCP, data_marker = LoadData(path_participant_data)