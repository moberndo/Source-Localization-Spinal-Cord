""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Calcualtion for Dipole Positions
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORTS -----
from math import pi
from ngsolve import Parameter
from numpy import linspace, array, append

# ----- FUNCTIONS -----
def DipolePositions(params_geom, num_dipoles):
    # dipoles: [x, y, z, alpha, beta, Parameter(p)]

    rad_dipole = params_geom[0][1] * 0.4

    dipoles_per_section = num_dipoles // 7
    num_dipoles = dipoles_per_section * 7

    #z_vals = linspace(start_dipole, end_dipole, num_dipoles)
    z_vals = array([])

    num_vert = params_geom[1][2]
    len_spine = params_geom[0][0]
    height_vert = params_geom[1][3]
    #height_thora = params_geom[1][5]
    height_disc_vert = params_geom[2][0]

    for i in range(num_vert):
        start_i = len_spine - (i+1) * (height_vert + height_disc_vert) + height_vert * 0.3
        end_i = start_i - height_vert * 0.3

        z_vals_i = linspace(start_i, end_i, dipoles_per_section)
        z_vals = append(z_vals, z_vals_i)


    #y_vals = [-rad_dipole, rad_dipole]
    y_vals = [rad_dipole]
    alpha_vals = [0, pi/2, 0]
    beta_vals = [0, 0, pi/2]

    dipoles = []
    for y in y_vals:
        for z in z_vals:
            dip = []
            for angle in range(3):
                dip.append([0, y, z, alpha_vals[angle], beta_vals[angle], Parameter(1)]) 
            dipoles.append(dip)

    print('Dipole positions calculated!')
    return dipoles, num_dipoles, z_vals
