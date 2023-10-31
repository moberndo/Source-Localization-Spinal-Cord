from numpy import zeros, save, identity, ones
from numpy import load
from numpy.linalg import pinv
from .fem_functions import CalcGFU


def LeadField(num_participant, main_path, leadfield_params):
    folder_participant = 'participant_' + num_participant
    leadfield_name = 'leadfield_NV' + str(len(leadfield_params[0])) + '_NE' + str(len(leadfield_params[1])) + '.npy'
    path_LeadField = main_path / 'data' / folder_participant / leadfield_name

    if path_LeadField.exists():
        K = load(path_LeadField)
        print('Existing Lead field found and loaded.')
    else:
        K = CalcLeadField(leadfield_params[0], leadfield_params[1], leadfield_params[2], path_LeadField)
        print('Lead field calculated!')

    return K


def CalcLeadField(dipoles, electrodes, fem_params, path):
    # The result is saved in the Lead Field K
    # shape(K) = N_E x N_V x 3
    #   - N_E : Number of electrodes
    #   - N_V : Number of dipoles
    #   - 3rd dimension in 3 because of the possible orientations of the dipole
    
    f = fem_params[0]
    maxh_ = fem_params[1]
    a = fem_params[2]
    c = fem_params[3]
    fes = fem_params[4]
    mesh = fem_params[5]
    # save NV NE to leadfiled name!!!!!!!!!!!
    N_V = len(dipoles)
    N_E = len(electrodes)
    K = zeros(shape = (N_E, N_V, 3))
    
    for idx_dip, dip in enumerate(dipoles):
        idx_dim = -1
        print(f'Dipole {idx_dip + 1} / {len(dipoles)}')
        for i in range(3):
            idx_dim += i
            dipole = dip[i]
            gfu, _ , _ = CalcGFU([dipole],f,maxh_, a, c, fes)
            for idx_elec, elec in enumerate(electrodes):
                elec_mesh = mesh(elec[0], elec[1], elec[2])
                gfu_res = gfu(elec_mesh)
                K[idx_elec, idx_dip, idx_dim] = gfu_res
    
    save(path, K)
    return K

def CalcCenteringMatrix(n):
    I = identity(n)
    vec_ones = ones((n,1))
    H = I - (vec_ones @ vec_ones.T) / (vec_ones.T @ vec_ones)
    return H

def CalcInvMatrix(K, H, alpha):
    p_inv = H @ K @ K.T @ H + alpha * H
    T = K.T @ H @ pinv(p_inv)
    return T






