""" ###################################################################
    Title: Solving the source localization problem of spinal cord 
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Evaluate Solution
    Author: Markus E. Oberndorfer
################################################################### """

# ----- IMPORT PACKAGES -----
from numpy.linalg import norm
from numpy import argmin

# ----- FUNCTIONS -----
def CalcFunctional(data_SCP, K, J, alpha):
    """
    Calcualtes the Error-Functional for a given solution vector and 
    spinal cord potential recordings (data_SCP).

    Input:
        - data_SCP: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data. [Numpy array]
        - K: Leadfield matrix, that stores the electrode potentials for every
        possible dipole setting. [Numpy array]
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]
        - alpha: Regularitation parameter for error funcational calculation

    Output:
        F: Error-Functional. [float]
    """
    F = norm(data_SCP - K @ J) + alpha * norm(J)
    return F

def GetBestParams(vals, alpha_vals, num_dip_vals):
    """
    For a given array of error functionals, the lowest value is searched in
    order to return the respective hyperparameters alpha and 
    num_dipoles_per_section.

    Input: 
        - vals: Numpy array that contains all the error functional values.
        [Numpy array]
        - alpha_vals: list of all alpha values that were used in the error
        functional calculations. [list, floats]
        - num_dip_vals: list of all num_dipoles_per_section values that
        were used in the error functional calculations. [list, floats]

    Output:
         - best_alpha_val: best alpha value to minimize error functional.
         [float]
         - best_num_dip_val: best num_dipoles_per_section value to
         minimize error functional. [int]
    """
    n_cols = len(num_dip_vals)

    best_idx = argmin(vals)
    row = best_idx // n_cols
    col = best_idx-n_cols*row

    return alpha_vals[row], num_dip_vals[col]

if __name__ == '__main__':
    ...
