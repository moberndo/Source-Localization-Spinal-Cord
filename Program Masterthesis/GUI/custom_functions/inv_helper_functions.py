""" #######################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Functions for FEM calculations
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from numpy import zeros, save, identity, ones, linspace
from numpy import load, reshape
from numpy.linalg import pinv
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ----- IMPORT CUSTOM FUNCTIONS -----
from .fem_functions import CalcGFU
from .dipole_positions import CalcDipoleFromLeadField

# ----- FUNCTIONS -----
def LeadField(mesh_folder, leadfield_params):
    """
    Function that eiter loads an exisitng leadfield or calls
    CalcLeadField() to calculate a new one.

    Input:
        - mesh_folder: Path to the mesh folder. [pathlib object]
        - main_path: Contains the path of the main file [pathlib object]
        - leadfield_params: List of parameters that are necessary for 
        LeadField calculations. [list]
    
    Output:
        - K: Leadfield matrix, that stores the electrode potentials for every
        possible dipole setting. [Numpy array]
    """
    leadfield_name = 'leadfield_NV' + str(len(leadfield_params[0])) \
                   + '_NE' + str(len(leadfield_params[1])) + '.npy'
    path_LeadField = mesh_folder / leadfield_name

    if path_LeadField.exists():
        K = load(path_LeadField)
        print('Existing Lead field found and loaded.')
    else:
        print('No Lead field found. Calculating now ...')
        K = CalcLeadField(leadfield_params[0], leadfield_params[1],
                          leadfield_params[2], path_LeadField)
        print('Lead field calculated!')

    return K

def CalcLeadField(dipoles, electrodes, fem_params, path, save_ = True):
    """
    If no pre-calculated leadfield is available, this function
    is called to compute a new one for the given parameters.

    Input:
        - dipoles: list of all possible dipoles that should be tested. Each
        list entry contains the location, direction and magnitude of one
        dipole. [list; list]
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode. [list, Numpy array]
        - fem_params: Parameters that are needed for FE calculations. [list]
        - path: Path at which the leadfield will be saved. [Pathlib object]

    Output:
        - K: Leadfield matrix, that stores the electrode potentials for every
        possible dipole setting. [Numpy array]
    """
    # The result is saved in the Lead Field K
    # shape(K) = N_E x (3 x N_V)
    #   - N_E : Number of electrodes
    #   - N_V : Number of dipoles
    
    # FEM params
    f = fem_params[0]
    maxh_ = fem_params[1]
    a = fem_params[2]
    c = fem_params[3]
    fes = fem_params[4]
    mesh = fem_params[5]

    # Size of the leadfield matrix
    N_V = len(dipoles)
    N_E = len(electrodes)
    K = zeros(shape = (N_E, N_V))

    # Computes the leadfield matrix in a  for loop
    for idx_dip, dip in enumerate(dipoles):
        print(f'Dipole {idx_dip + 1} / {len(dipoles)}')
        gfu, _ , _ = CalcGFU([dip],f,maxh_, a, c, fes)
        for idx_elec, elec in enumerate(electrodes):
            elec_mesh = mesh(elec[0], elec[1], elec[2])
            gfu_res = gfu(elec_mesh)
            K[idx_elec, idx_dip] = gfu_res
    # Save the leadfield matrix to a given path
    if save_:
        save(path, K)

    return K

def CalcCenteringMatrix(n):
    """
    Computes the centering matrix that is needed for the sLORETA algorithm

    Input: 
        - n: size of the centering matrix. [int]

    Output:
        - H: Centering matrix of size (n x n). [Numpy array]
    """
    I = identity(n)
    vec_ones = ones((n,1))
    H = I - (vec_ones @ vec_ones.T) / (vec_ones.T @ vec_ones)
    return H

def CalcInvMatrix(K, H, alpha):
    """
    Computes the inverse of the leadfield matrix according to the
    sLORETA algorithm.

    Input:
        - K: Leadfield matrix, that stores the electrode potentials for every
        possible dipole setting. [Numpy array]
        - H: Centering matrix of size (n x n). [Numpy array]
    
    Output:
        - T: Inverse of the leadfield matrix. [Numpy array]
    """
    p_inv = H @ K @ K.T @ H + alpha * H
    T = K.T @ H @ pinv(p_inv)
    return T

def CalcSimulatedPotentials(J, electrodes, fem_params, params_geom,
                            num_dipoles):
    """
    For a given solution vector, this function computes the potential that
    would occur at the electrodes.

    Input: 
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode. [list, Numpy array]
        - fem_params: Parameters that are needed for FE calculations. [list]
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
        - num_dipoles: Total number of dipoles. [int]

    Output:
        - data_sim: Contains the potentials at the electrods caused by
        the given solution vector J. [Numpy array]  
    """

    f = fem_params[0]
    maxh_ = fem_params[1]
    a = fem_params[2]
    c = fem_params[3]
    fes = fem_params[4]
    mesh = fem_params[5]

    num_dipoles_per_section = int(num_dipoles / 7 / 3)

    dipoles, _, _ = CalcDipoleFromLeadField(params_geom, J,
                                            num_dipoles_per_section)
    J = reshape(J, newshape=(1, J.shape[0]))
    len_signal = J.shape[0]
    N_E = len(electrodes)
    data_sim = zeros(shape = (N_E, len_signal))
    N_V = len(dipoles[0])
    
    for idx_t in range(len_signal):
        print(f'Signal point {idx_t + 1} / {len_signal}')
        all_dipoles = []
        for idx_dip in range(N_V):
            for i in range(3):
                #if i == 2:
                all_dipoles.append(dipoles[idx_t][idx_dip][i])
        gfu, _ , _ = CalcGFU(all_dipoles,f,maxh_, a, c, fes)

        for idx_elec, elec in enumerate(electrodes):
            elec_mesh = mesh(elec[0], elec[1], elec[2])
            gfu_res = gfu(elec_mesh)
            data_sim[idx_elec, idx_t] = gfu_res
    
    return data_sim