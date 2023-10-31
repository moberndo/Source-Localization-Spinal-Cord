
from math import pi
from ngsolve import Parameter
from numpy import arange, sin, ones, save, load, amax
from numpy.random import rand
from .fem_functions import CalcGFU, FEM_Params




def CreateArtificialSignal(dipole_loc, electrodes, mesh, mat_conduc, maxh_, params_geom, signal_length, main_path,
                           SNR = None, dir = 'x', side = 'left', fs = 512):
    # creates an artificial potential distribution caused by 
    #   - case 1: 1 dipole (one location)
    #   - case 2: >1 neighboring dipoles (one location)
    #   - case 3: >M dipoles (M locations)
    # the magnitude of the dipole(s) follows a sinus with max magnitude beeing 1 multiplied with some scalar 0.5 < s < 1.5
    
    
    def CreateDipole(x, y, z, alpha, beta, p):
        dip = [x,y,z,alpha, beta, Parameter(p)]
        return dip
    
    def CreateSignal(signal_length, fem_params, rad_dipole, y, electrodes, mesh, side_, dipole_loc,
                     alpha_beta, SNR = None):
        arti_sig = ones((16, signal_length))
        for idx, val in enumerate(y):
            print(f'Calculating Dipole {idx + 1} / {signal_length}')
            dip = [CreateDipole(rad_dipole * side_, 0, dipole_loc[0][0], alpha_beta[0], alpha_beta[1], val)]

            gfu, _ , _ = CalcGFU(dip,fem_params[0],fem_params[1], fem_params[2], fem_params[3], fem_params[4])

            for idx_elec, elec in enumerate(electrodes):
                elec_mesh = mesh(elec[0], elec[1], elec[2])
                gfu_res = gfu(elec_mesh)
                arti_sig[idx_elec,idx] = gfu_res

        if SNR != None:
            noise_level = amax(arti_sig) / SNR
            noise = rand(arti_sig.shape) * noise_level
            arti_sig += noise

        return arti_sig
    
    a, c, f, fes = FEM_Params(mesh, mat_conduc)
    fem_params = [f, maxh_ , a, c, fes, mesh]
    
    t = arange(signal_length) / fs
    f = 1 # frequency of the signal
    rad_dipole = params_geom[0][1] * 0.75

    if side == 'left':
        side_ = 1
    elif side == 'right':
        side_ = -1

    if dir == 'x':
        alpha_beta = [0, 0]
    elif dir == 'y':
        alpha_beta = [pi/2, 0]
    elif dir == 'z':
        alpha_beta = [0, pi/2]

    mesh = fem_params[5]

    # read dipole location
    

    if len(dipole_loc) == 1:
        if len(dipole_loc[0]) == 1:
            # case 1
            name_file = 'artificial_sig_' + str(len(dipole_loc)) + '_' + str(dipole_loc[0][0]) + '_' + dir + '_' + side + '_' + str(fs) + '_' + str(signal_length)
            if SNR == None:
                name_file = name_file + '.npy'
            else: 
                name_file = name_file + '_SNR' + str(SNR) + '.npy'

            path_arti_sig = main_path / 'data' / 'artificial_data' / name_file

            if path_arti_sig.exists():
                arti_sig = load(path_arti_sig)
                return arti_sig
            else:
                y = sin(2* pi * f * t)
                arti_sig = CreateSignal(signal_length, fem_params, rad_dipole, y, electrodes, mesh,
                                        side_, dipole_loc, alpha_beta, SNR = None)
                save(path_arti_sig, arti_sig)
                return arti_sig

        elif len(dipole_loc[0]) > 1:
            # case 2
            ...
    elif len(dipole_loc) > 1:
        # case 3
        ...

    
    return arti_sig

if __name__ == '__main__':
    #arti_sig = CreateArtificialSignal(dipole_loc, electrodes, fem_params, params_geom, signal_length)
    ...