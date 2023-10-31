""" #######################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Functions for Inverse Probelm Calculations
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from numpy import ones, reshape

# ----- IMPORT CUSTOM FUNCTIONS -----
from .inv_helper_functions import CalcCenteringMatrix, CalcInvMatrix
from .inv_helper_functions import LeadField
from .fem_functions import FEM_Params

# ----- FUNCTIONS -----
def sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles, electrodes,
            mesh_folder, main_path, alpha):
    """
    Calcualte the inverse problem for a given geometry and spinal cord
    potential recordings and several other parameters. (see below)

    Input:
        - data_SCP: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data. [Numpy array]
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - mat_conduc: Dictionary that stores the conductivity values for all
         materials in the model. [dict]
        - maxh_: Is the maximum size of one mesh element. [float]
        - dipoles: list of all possible dipoles that should be tested. Each
        list entry contains the location, direction and magnitude of one
        dipole. [list; list]
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode. [list, Numpy array]
        - mesh_folder: Path to the mesh folder. [pathlib object]
        - main_path: Contains the path of the main file [pathlib object]
        - alpha: alpha is the regularization parameter. [float]

    Output: 
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]
    """
    # Calculate constants
    N_E = len(electrodes)
    N_D = len(dipoles)

    signal_length = data_SCP.shape[1]

    # Calculate parameters
    a, c, f, fes = FEM_Params(mesh, mat_conduc)
    fem_params = [f, maxh_ , a, c, fes, mesh]
    leadfield_params = [dipoles, electrodes, fem_params]

    # Calculate Matrices
    K = LeadField(mesh_folder, leadfield_params)
    H = CalcCenteringMatrix(N_E)

    J = ones(shape = (signal_length,3 * N_D))

    T = CalcInvMatrix(K, H, alpha = alpha)

    # Solve inverse problem
    J = ones(shape = (signal_length,N_D))
    for j in range(signal_length):
        J_ = T @ data_SCP[:,j]
        J[j,:] = reshape(J_, newshape = (1, N_D))

    return J