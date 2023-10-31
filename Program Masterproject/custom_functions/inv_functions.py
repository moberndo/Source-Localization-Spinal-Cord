# imports
from .inv_helper_functions import CalcCenteringMatrix, CalcInvMatrix, LeadField
from .fem_functions import FEM_Params
from numpy import ones, reshape

def sLORETA(data_SCP, mesh, mat_conduc, maxh_, dipoles, electrodes, num_participant, main_path):

    # Calculate constants
    N_E = len(electrodes)
    N_D = len(dipoles)
    signal_length = data_SCP.shape[1]

    # Calculate parameters
    a, c, f, fes = FEM_Params(mesh, mat_conduc)
    fem_params = [f, maxh_ , a, c, fes, mesh]
    leadfield_params = [dipoles, electrodes, fem_params]

    # Calculate Matrices
    K = LeadField(num_participant, main_path, leadfield_params)
    H = CalcCenteringMatrix(N_E)

    J = ones(shape = (signal_length,N_D,3))
    for i in range(3):
        T_i = CalcInvMatrix(K[:,:,i], H, alpha = 0.005)

        # Solve inverse problem
        J_i = ones(shape = (signal_length,N_D))
        for j in range(signal_length):
            J_ij = T_i @ data_SCP[:,j]
            J_i[j,:] = reshape(J_ij, newshape = (1,N_D))
        J[:,:,i] = J_i
    return J