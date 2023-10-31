""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Initilization file for custom functions
    Author: Markus E. Oberndorfer
####################################################################### """


from .fem_functions import FEM_Params, CalcGFU
from .mesh_functions import create_mesh_1, create_mesh_2, create_mesh_3, CreateParamsGeom
from .inv_helper_functions import CalcLeadField, CalcCenteringMatrix, LeadField, CalcInvMatrix
from .load_data import LoadData
from .process_data import ProcessData
from .dipole_positions import DipolePositions
from .electrode_positions import ElectrodePositions
from .model_selection import ModelSelection
from .inv_functions import sLORETA
from .visualize_solution import DrawBackground, DrawSolution, CreateAnimation
from .create_artificial_signal import CreateArtificialSignal