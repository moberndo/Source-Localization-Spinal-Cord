""" ###################################################################
    Title: Solving the source localization problem of spinal cord 
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Calculate Hyperparameters
    Author: Markus E. Oberndorfer
################################################################### """

# ----- IMPORT PACKAGES -----
from numpy import save, load, zeros, meshgrid, mean
from os.path import isfile
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import cm

# ----- IMPORT CUSTOM FUNCTIONS -----
from .inv_functions import sLORETA
from .inv_helper_functions import LeadField
from .visualize_solution import GetWaveformPoints
from .fem_functions import FEM_Params
from .dipole_positions import DipolesLeadfield
from .evaluate import CalcFunctional, GetBestParams

# ----- FUNCTIONS -----
def HyperParameterSearch(wp_type, alpha_hp, dipoles_hp, electrodes,
                         user_ID, main_path, data_SCP, mesh,
                         mat_conduc, maxh_, fs):
    """
    For a given waveform-point (wp)-type and a list of hyperparameters, the
    corresponding error functional F is calcualted.

    Input:
        - wp_type: What waveform points should be calculated. Possible types
        are P1, N1 and P2. [str]
        - alpha_hp: list of all possible alpha values that should be tested.
        alpha is the regularization parameter. [list; float]
        - dipoles_hp: list of all possible dipoles that should be tested. Each
        list entry contains the location, direction and magnitude of one
        dipole. [list; list]
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode [list, Numpy array]
        - user_ID: Is a short ID to points to the folder structure of
        interest. [str]
        - main_path: Contains the path of the main file [pathlib object]
        - data_SCP: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data. [Numpy array]
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - mat_conduc: Dictionary that stores the conductivity values for all
         materials in the model. [dict]
        - maxh_: Is the maximum size of one mesh element. [float]
        - fs: Sampling rate during the measurment. [int]

    Output:
        - F_vals: Is an array of shape (LengthAlphaList x LengthDipolesList)
        that contains every error functional result. [Numpy array]
    """
    f_vals_path = 'F_vals_' + wp_type + '.npy'
    user_folder = 'user_' + user_ID
    f_vals_path = main_path / 'data' / user_folder / f_vals_path
    if not isfile(f_vals_path):
        F_vals = zeros((len(alpha_hp), len(dipoles_hp)))
        for idx_alpha, alpha in enumerate(alpha_hp):
            for idx_dipoles, dipoles in enumerate(dipoles_hp):
                J = sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles,
                            electrodes, user_ID, main_path, alpha)
                P1, N1, P2 = GetWaveformPoints(J, data_SCP, fs)
                if wp_type == 'P1':
                    J_wp = J[P1,:]
                    data_SCP_wp = data_SCP[:,P1]
                elif wp_type == 'N1':
                    J_wp = J[N1,:]
                    data_SCP_wp = data_SCP[:,N1]
                elif wp_type == 'P2':
                    J_wp = J[P2,:]
                    data_SCP_wp = data_SCP[:,P2]
                else:
                    print('Waveform point is not valid.')
                a, c, f, fes = FEM_Params(mesh, mat_conduc)
                fem_params = [f, maxh_, a, c, fes, mesh]

                leadfield_params = [dipoles, electrodes, fem_params]
                K = LeadField(user_ID, main_path, leadfield_params)
                func = (CalcFunctional(data_SCP_wp, K, J_wp , alpha))
                F_vals[idx_alpha, idx_dipoles] = func

        save(f_vals_path, F_vals)
    else:
        F_vals = load(f_vals_path)
    return F_vals

def CalcHyperparams(params_geom, num_dips_hp, alpha_hp, electrodes,
                    user_ID, main_path, data_SCP, mesh, mat_conduc,
                    maxh_, fs, plot_data = False):
    """
    For a given waveform-point (wp)-type and a list of hyperparameters, the
    corresponding error functional F is calcualted, from which the best Hyper-
    parameters are chosen and returned. If plot_data is set to True, all
    values for F are returned, not just the hyperparameters

    Input:
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
        - num_dips_hp: list of number of dipoles per intervertebral
          section that are assumed a priori. [list; int]
        - alpha_hp: list of all possible alpha values that should be tested.
        alpha is the regularization parameter. [list; float]
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode [list, Numpy array]
        - user_ID: Is a short ID that points to the folder structure of
        interest. [str]
        - main_path: Contains the path of the main file [pathlib object]
        - data_SCP: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data. [Numpy array]
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - mat_conduc: Dictionary that stores the conductivity values for all
         materials in the model. [dict]
        - maxh_: Is the maximum size of one mesh element. [float]
        - fs: Sampling rate during the measurment. [int]
        - plot_data: Boolean variable that, if set to true, returns all error-
        functional values. [bool]

    Output:
        - list_1: List with two entries. First entry is the best alpha value,
        second entry is the best num_dipoles_per_section value.
        [list; float & int]
        - list_2 (optional): If plot_data is True, returns F_vals for all
        waveform point types. For more information on F_vals type
        help(HyperParameterSearch). [list; Numpy array]
    """
    
    # calc hp params
    dipoles_hp = [DipolesLeadfield(params_geom, num_dips)[0]
                  for num_dips in num_dips_hp]

    F_vals_P1 = HyperParameterSearch('P1', alpha_hp, dipoles_hp, electrodes,
                                     user_ID, main_path, data_SCP,
                                     mesh, mat_conduc, maxh_, fs)
    
    best_param_idx_P1 = GetBestParams(F_vals_P1, alpha_hp, num_dips_hp)
    alpha_P1 = best_param_idx_P1[0]
    num_dip_P1 = best_param_idx_P1[1]

    F_vals_N1 = HyperParameterSearch('N1', alpha_hp, dipoles_hp, electrodes,
                                     user_ID, main_path, data_SCP,
                                     mesh, mat_conduc, maxh_, fs)
    
    best_param_idx_N1 = GetBestParams(F_vals_N1, alpha_hp, num_dips_hp)
    alpha_N1 = best_param_idx_N1[0]
    num_dip_N1 = best_param_idx_N1[1]

    F_vals_P2 = HyperParameterSearch('P2', alpha_hp, dipoles_hp, electrodes,
                                     user_ID, main_path, data_SCP,
                                     mesh, mat_conduc, maxh_, fs)

    best_param_idx_P2 = GetBestParams(F_vals_P2, alpha_hp, num_dips_hp)
    alpha_P2 = best_param_idx_P2[0]
    num_dip_P2 = best_param_idx_P2[1]

    alpha = mean([alpha_P1, alpha_N1, alpha_P2])
    num_dip = round(mean([num_dip_P1, num_dip_N1, num_dip_P2]))

    if plot_data:
        return [alpha, num_dip], [F_vals_P1, F_vals_N1, F_vals_P2]
    else:
        return [alpha, num_dip]


def VisualizeHyperParams(F_vals, num_dips_hp, alpha_hp, wp_type, save = False,
                          main_path = None, user_ID = None):
    
    """
    To visualize the impact of the hyperparameters on the error functional,
    this function created a 3D surface plot. If save is set to True, the 
    variables main_path and user_ID have to be set to create the image path.

    Input:
        - F_vals: Is an array of shape (LengthAlphaList x LengthDipolesList)
        that contains every error functional result. [Numpy array]
        - num_dips_hp: list of number of dipoles per intervertebral
          section that are assumed a priori. [list; int]
        - alpha_hp: list of all possible alpha values that should be tested.
        alpha is the regularization parameter. [list; float]
        - wp_type: What waveform points should be calculated. Possible types
        are P1, N1 and P2. [str]
        - save: If true, the image is saved. [bool]
        - main_path: Contains the path of the main file [pathlib object]
        - user_ID: Is a short ID to points to the folder structure of
        interest. [str]

    Output:
        - No output.
    """
    X = num_dips_hp
    Y = alpha_hp
    X, Y = meshgrid(X, Y)
    Z = F_vals

    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams['axes.titlepad'] = -6
    rcParams["xtick.major.size"] = 0
    

    _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.view_init(elev=30, azim=125, roll=0)
    
    ax.set_ylabel(r'$\mathtt{alpha}$')
    ax.set_xlabel(r'$\mathtt{num\_dipoles\_per\_section}$')
    ax.set_zlabel('error functional F', rotation=90)
    ax.set_facecolor('white')
    #ax.w_xaxis.pane.fill = False
    #ax.w_yaxis.pane.fill = False
    #ax.w_zaxis.pane.fill = False
    ax.set_yticks([0.005, 0.01, 0.05, 0.1, 0.2, 0.35, 0.5],
                  ['', '0.01', '', '0.1', '0.2', '0.35', '0.5'])
    ax.set_xticks([1,2,3,4,5], ['1','2','3','4','5'])

    ax.plot_surface(X, Y, Z, cmap=cm.Reds)
    title_ = 'hyperparameter search for: ' + wp_type
    ax.set_title(title_)

    if save:
        user_folder = 'user_' + user_ID
        file_name = 'hp_search_' + wp_type + '.png'
        file_path = main_path / 'data' / user_folder / 'images' / file_name
        plt.savefig(file_path)

if __name__ == '__main__':
    ...