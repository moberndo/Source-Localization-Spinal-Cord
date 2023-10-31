""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate
      the voltage distribution caused by a dipole
    Subtitle: Initilization file for custom functions
    Author: Markus E. Oberndorfer
####################################################################### """


from .fem_functions import FEM_Params, CalcGFU
from .mesh_functions import CreateMesh, CreateParamsGeom
from .inv_helper_functions import CalcLeadField, CalcCenteringMatrix, \
                                  LeadField, CalcInvMatrix, \
                                  CalcSimulatedPotentials
from .load_data import LoadData
from .process_data import ProcessData, ChannelPlotWP_1, ChannelPlotWP_2
from .dipole_positions import DipolesLeadfield, CalcDipoleFromLeadField
from .electrode_positions import ElectrodePositions
from .model_selection import ModelSelection
from .inv_functions import sLORETA
from .visualize_solution import ShowSolution, GetWaveformPoints
from .evaluate import CalcFunctional
from .calc_hyperparams import CalcHyperparams, VisualizeHyperParams, \
                              HyperParameterSearch
