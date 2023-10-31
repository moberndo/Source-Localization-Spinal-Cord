"""###########################################################################
        TITLE:  Graphical User Interface for Source Localization Algorithm
        AUTHOR: Markus E. Oberndorfer
###########################################################################"""

# ---- IMPORT PACKAGES ----
from PyQt5.QtWidgets import (QMainWindow, QApplication, QLabel, QLineEdit,
                             QPushButton, QRadioButton, QComboBox,
                             QFileDialog, QMessageBox)
from PyQt5 import uic
import sys
import os
from os import mkdir
from os.path import isdir
from pathlib import Path

# ---- IMPORT CUSTOM FUNCTIONS ----

from custom_functions import ModelSelection, LoadData, ProcessData
from custom_functions import DipolesLeadfield, ElectrodePositions
from custom_functions import sLORETA, ShowSolution, GetWaveformPoints

# ---- IMPORT PARAMETERS ----
from params import params_elec, mat_conduc

# ---- MAIN ----

class UI(QMainWindow):
    def __init__(self):
        super(UI, self).__init__()
        self.vtk_widget = None

        # Load UI file
        uic.loadUi('interface.ui', self)

        # Define Widgets in 'Select EMG Data'
        self.rdb_SingleParticipant = self.findChild(QRadioButton, 'rdbtn_SinglePart')
        self.rdb_SingleParticipant.setChecked(True)
        self.rdb_MultiParticipant = self.findChild(QRadioButton, 'rdbtn_MultiPart')
        self.rdb_MultiParticipant.setChecked(False)
        self.btn_LoadData = self.findChild(QPushButton, 'btn_LoadData')
        self.lbl_DataPath = self.findChild(QLabel, 'lbl_DataPath')

        # Define Widgets in 'Choose Geometry'
        self.rdb_AltGeomParams = self.findChild(QRadioButton, 'rdbtn_AltGeom')
        self.rdb_AltGeomParams.setChecked(True)
        self.rdb_NewGeomParams = self.findChild(QRadioButton, 'rdbtn_GeomMeasur')
        self.rdb_NewGeomParams.setChecked(False)
        self.ledit_InionC7 = self.findChild(QLineEdit, 'ledit_InionC7')
        self.ledit_C7L5 = self.findChild(QLineEdit, 'ledit_C7L5')
        self.ledit_NeckCircum = self.findChild(QLineEdit, 'ledit_NeckCircum')

        # Define Widgets in 'Methods'
        self.rdb_sLORETA = self.findChild(QRadioButton, 'rdbtn_sLORETA')
        self.rdb_sLORETA.setChecked(True)
        self.cbox_NumElec = self.findChild(QComboBox, 'cBox_NumElec')
        self.cbox_NumDips = self.findChild(QComboBox, 'cbox_NumDipoles')

        # Define Widgets in 'Start'
        self.btn_StartCalc = self.findChild(QPushButton, 'btn_StartCalc')

        # Define Widghts in 'Results'
        self.lbl_Feedback = self.findChild(QLabel, 'lbl_Feedback')

        # Load the data
        self.btn_LoadData.clicked.connect(self.load_data)

        # Start the calculation
        self.btn_StartCalc.clicked.connect(self.start_calc)

        # Show the Application
        self.show()

    def show_messagebox(self, title, message, style):
        msg = QMessageBox()
        msg.setWindowTitle(title)
        msg.setText(message)
        if style == 'i':
            msg.setIcon(QMessageBox.Information)
        elif style == 'q':
            msg.setIcon(QMessageBox.Question)
        elif style == 'w':
            msg.setIcon(QMessageBox.Warning)

        x = msg.exec_()

    def load_data(self):
        # Check selection:
        if not self.rdb_SingleParticipant.isChecked() and not self.rdb_MultiParticipant.isChecked():
            title = 'Information'
            message = 'You need to select an option. Select "Single participant" or "Multiple participants".'
            self.show_messagebox(title, message, 'i')
        else:
            folder_path = QFileDialog.getExistingDirectory(self, 'Select Folder')
            files = os.listdir(str(folder_path))
            cancel = False

            if self.rdb_SingleParticipant.isChecked():
                for file in files:
                    _, ext = os.path.splitext(file)
                    if ext != '.mat':
                        title = 'Warning'
                        message = 'Your selected folder should only contain .mat files.'
                        self.show_messagebox(title, message, 'w')
                        cancel = True
                        break
                if not cancel:
                    self.lbl_DataPath.setText(str(folder_path))

            elif self.rdb_MultiParticipant.isChecked():
                for file in files:
                    _, ext = os.path.splitext(file)
                    if ext != '':
                        title = 'Warning'
                        message = 'Your selected folder should only contain other folders. All subfolders should only contain .mat files.'
                        self.show_messagebox(title, message, 'w')
                        cancel = True
                        break
                if not cancel:
                    self.lbl_DataPath.setText(str(folder_path))
        
    def check_settings(self):
            def isnumber(string1):
                if string1.isnumeric():
                    return True
                else:
                    try:
                        float(string1)
                        return True
                    except:
                        return False

            if os.path.isdir(str(self.lbl_DataPath.text())):
                if not self.rdb_NewGeomParams.isChecked() and not self.rdb_AltGeomParams.isChecked():
                    title = 'Information'
                    message = 'You need to select a geometry-setting. Select "Alternative Geometry Parameters" or "New Geometry Parameter".'
                    self.show_messagebox(title, message, 'i')
                    return False
                
                if self.rdb_NewGeomParams.isChecked():
                    if not (self.ledit_InionC7.text() != '' and isnumber(self.ledit_InionC7.text())):
                        title = 'Information'
                        message = 'You need to set all parameters.'
                        self.show_messagebox(title, message, 'i')
                        return False
                    if not (self.ledit_C7L5.text() != '' and isnumber(self.ledit_C7L5.text())):
                        title = 'Information'
                        message = 'You need to set all parameters.'
                        self.show_messagebox(title, message, 'i')
                        return False
                    if not (self.ledit_NeckCircum.text() != '' and isnumber(self.ledit_NeckCircum.text())):
                        title = 'Information'
                        message = 'You need to set all parameters.'
                        self.show_messagebox(title, message, 'i')
                        return False
                    
                if not self.rdb_sLORETA.isChecked():
                        title = 'Information'
                        message = 'You need to select a method.'
                        self.show_messagebox(title, message, 'i')
                        return False
                
            return True

    def start_calc(self):
        
        if self.check_settings():
            # Calculate Solution
            J, params_geom, z_vals, num_dipoles_per_section, data_SCP, fs = self.main()
            # Visualize Solution
            P1, N1, P2 = GetWaveformPoints(J, data_SCP, fs)
            J_WP = [J[P1,:], J[N1,:], J[P2,:]]
            ShowSolution(J_WP, params_geom, z_vals, num_dipoles_per_section)

        else:
            title = 'Warning'
            message = 'Some settings are not valid.'
            self.show_messagebox(title, message, 'w')

    def main(self):   
        main_path = Path(__file__).parent
        print(str(main_path))
        if self.rdb_AltGeomParams.isChecked():
            alt_path = main_path / 'data' / 'measurment_geom_alt.txt'
            if alt_path.exists():
                geom = open(alt_path, 'r')
                measurments_geom = geom.read()
                geom.close()
                meas_geom_list = measurments_geom.split(';')
                Inion_C7 = meas_geom_list[0]
                C7_L5 = meas_geom_list[1]
                ledit_NeckCircum = meas_geom_list[2]
            else:
                exit('File missing! Alternative geometry measure not available.')
        else:
            Inion_C7 = str(self.ledit_InionC7.text())
            C7_L5 = str(self.ledit_C7L5.text())
            ledit_NeckCircum = str(self.ledit_NeckCircum.text())
            measurments_geom = Inion_C7 + ';' + C7_L5 + ';' + ledit_NeckCircum
        print(measurments_geom)
    
        num_electrodes = int(self.cbox_NumElec.currentText())
        num_dipoles_per_section = int(self.cbox_NumDips.currentText())
        fs = 512
        alpha = 0.05
        data_path = self.lbl_DataPath.text()
        #------------------------------------#

        ################## MAIN ###################
        # Create folder structure
        mesh_folder = 'mesh_' + ledit_NeckCircum + '_' + Inion_C7 + '_' + C7_L5
        mesh_folder = main_path / 'data' / mesh_folder
        if not isdir(mesh_folder):
            mkdir(mesh_folder)

        #########################################
        '''' ======== LOAD ALL FILES ======== '''
        #########################################

        # Select the model or load existing model
        mesh, params_geom, loc_C7, maxh_ = ModelSelection(measurments_geom,
                                                          main_path,
                                                          mesh_folder)
        print('Model sucessfully loaded.')
        self.lbl_Feedback.setText('Model sucessfully loaded.')

        # Load the data and process it
        if self.rdb_SingleParticipant.isChecked():
            type = 'single'
        else:
            type = 'multiple'
        _ , data_SCP, data_marker, num_participants = LoadData(data_path, type)
        print('Data loaded.')
        self.lbl_Feedback.setText('Data sucessfully loaded.')
        data_SCP, _ = ProcessData(data_SCP, data_marker, num_participants,
                                  fs, main_path)
        print('Data processed.')
        self.lbl_Feedback.setText('Data sucessfully processed.')
        

        # Define the position of the electrodes (depending on C7)
        electrodes, N_E = ElectrodePositions(loc_C7, params_elec, params_geom)
        print('Electrode positions calculated.')

        #########################################
        '''' ==== SOLVE INVERSE PROBLEM ==== '''
        #########################################

        # Create dipoles that are necessary for the Leadfield-Calculation
        dipoles, z_vals = DipolesLeadfield(params_geom, num_dipoles_per_section)
        print('Dipoles calculated.')
        # Solve the inverser problem
        J = sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles, electrodes, 
                    mesh_folder, main_path, alpha)

        print('Inverse Problem calculated.')
        self.lbl_Feedback.setText('Inverse Problem calculated.')
        return J, params_geom, z_vals, num_dipoles_per_section, data_SCP, fs

        #####################################
        '''' ==== VISUALIZE SOLUTION ==== '''
        #####################################

        '''# Calculate the Waveform-Points, at which the solution is visualized
        P1, N1, P2 = GetWaveformPoints(J, data_SCP, fs)
        J_WP = [J[P1,:], J[N1,:], J[P2,:]]

        # Create a 3D image of the solution at all time points P1, N1 and P2
        # view options:
        # frontal, lower_frontal, upper_frontal, upper_right, upper_left 
        ShowSolution(J_WP, params_geom, z_vals, num_dipoles_per_section,
                    view = 'frontal', save = True, path = images_folder)

        plt.show()'''


if __name__ == '__main__':
    app = QApplication(sys.argv)
    UIWindow = UI()
    app.exec_()