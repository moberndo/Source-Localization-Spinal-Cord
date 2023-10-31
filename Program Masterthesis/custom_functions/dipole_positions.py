""" ################################################################### 
    Title: Solving the source localization problem of spinal cord
           potentials for electric dipoles in the cervical spinal cord
    Subtitle: Calculation of Dipole Position
    Author: Markus E. Oberndorfer
################################################################### """

# ----- IMPORT PACKAGES -----
from math import pi
from ngsolve import Parameter
from numpy import linspace, array, append

# ----- FUNCTIONS -----
def DipolePositions(params_geom, num_dipoles_per_section):
    """
    DipolePosition takes the geometry parameters and the number of
    dipoles per section as input and outputs the respective 3d 
    coordinates for this geometry

    Input:
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
        - num_dipoles_per_section: number of dipoles per intervertebral
          section that are assumed a priori. [int]

    Output:
        - y_vals: y values for all dipoles [list; float]
        - z_vals: z values for all dipoles [list; float]
        - alpha_vals: alpha values for all dipoles. alpha is the angle
         in the transverse plane [list; float]
        - beta_vals: beta values for all dipoles. beta is the angle
        in the frontal plane [list; float]
        - num_dipoles: is the total number of dipoles [int]
    """
    rad_dipole = params_geom[0][1] * 0.25

    if num_dipoles_per_section < 1:
        raise ValueError('Number of dipoles per section has \
                          to be greater or equal than one.')
    
    num_dipoles = num_dipoles_per_section * 7
    z_vals = array([])

    num_vert = params_geom[1][2][0]
    len_spine = params_geom[0][0]
    height_vert = params_geom[1][2][5]
    height_disc_vert = params_geom[2][1]

    for i in range(num_vert):
        start_i = len_spine - (i+1) * (height_vert + height_disc_vert) \
                + height_vert * 0.3
        end_i = start_i - height_vert * 0.3

        z_vals_i = linspace(start_i, end_i, num_dipoles_per_section)
        z_vals = append(z_vals, z_vals_i)

    y_vals = [rad_dipole]
    
    alpha_vals = [0, pi/2, 0]
    beta_vals = [0, 0, pi/2]
    return y_vals, z_vals, alpha_vals, beta_vals, num_dipoles

def DipolesLeadfield(params_geom, num_dipoles_per_section):
    """
    DipolesLeadfield can be used to calculate the dipole
    parameters (location and direction) from a given for a
    given number of dipoles per intervertebral section.

    Input:
        - params_geom: software intern structure. See ReadMe.txt file
        for more information. [-]
        - num_dipoles_per_section: number of dipoles per intervertebral
        section that are assumed a priori. [int]

    Output:
        - dipoles: list of all dipoles with their respective location
        and orientation. [list; list]
        - z_vals: z values for all dipoles [list; int]
    """
    out_dipole_pos = DipolePositions(params_geom, num_dipoles_per_section)
    y_vals = out_dipole_pos[0]
    z_vals = out_dipole_pos[1]
    alpha_vals = out_dipole_pos[2]
    beta_vals = out_dipole_pos[3]
    print('Dipole positions calculated!')
    dipoles = []
    num_angles = len(alpha_vals)

    for y in y_vals:
        for z in z_vals:
            for angle in range(num_angles):
                dipoles.append([0, y, z, alpha_vals[angle],
                                beta_vals[angle], Parameter(1)]) 

    print('Dipoles calculated')
    return dipoles, z_vals

def CalcDipoleFromLeadField(params_geom, J, num_dipoles_per_section):
    """
    CalcDipoleFromLeadField can be used to calculate the dipole
    parameters (location and direction) from a given for a
    given number of dipoles per intervertebral section and solution
    vector, so that the dipoles have a specific magnitude.

    Input:
        - params_geom: software intern structure. See ReadMe.txt file
        for more information. [-]
        - J: Solution vector that contains the magnitude for every dipole.
        [list, float]
        - num_dipoles_per_section: number of dipoles per intervertebral
        section that are assumed a priori. [int]

    Output:
        - dipoles: list of all dipoles with their respective location
        and orientation. [list; list]
        - num_dipoles: Total number of dipoles. [int]
        - z_vals: z values for all dipoles [list; int]
    """
    J = J.reshape((-1,int(num_dipoles_per_section*7),3))
    len_signal = J.shape[0]
    num_dipoles_per_section = int(J.shape[1] / 7)
    num_angles = J.shape[2]

    y_vals, z_vals, alpha_vals, beta_vals, num_dipoles = DipolePositions(params_geom, num_dipoles_per_section)
    #print('Dipole positions calculated for Hyperparameter tuning!')
    dipoles = []
    y = y_vals[0]

    for i in range(len_signal):
        dipoles_t_i = []
        for idx_z, z in enumerate(z_vals):
            dip = []
            for idx_angle in range(num_angles):
                dip.append([0, y, z, alpha_vals[idx_angle], beta_vals[idx_angle], Parameter(J[i,idx_z, idx_angle])]) 
            dipoles_t_i.append(dip)
        dipoles.append(dipoles_t_i)
    print('Dipoles calculated from Lead Field!')
    return dipoles, num_dipoles, z_vals

if __name__ == '__main__':
    help(DipolePositions)
