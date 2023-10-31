""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Calcualtion for Electrode Positions
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORTS -----
from math import sqrt, sin, cos
from numpy import array

# ----- FUNCTIONS -----
def ElectrodePositions(electrode_pattern, loc_C7, params_elec, params_geom, show = False):
    
    dist_C7_Sp14 = params_elec[0]
    rad_electrode = params_elec[1]
    ver_dist_elec = params_elec[2]
    rad_neck = params_geom[3][0]
    loc_neck = params_geom[3][1]
    # pattern 1 is the electrode pattern from the paper 'wimmer, kostoglu and müller-putz'
    # pattern 2 is from the experimental setup of '...'
    if electrode_pattern == '16':
        ''' 
        Three columns on one side of the spine: 16 electrodes per side --> 32 eletrodes in total
        One 16 elec. array consists of:
            Left Col.: 5 electrodes
            Middle Col.: 6 electrodes
            Right Col.: 5 electrodes
        Vertical distance between electrodes (center to center) is 18mm
        '''

        # Correction for curved neck surface
        corr_x = lambda d, r: r - r * cos( d / r )
        corr_y = lambda d,r: r * sin( d / r)

        # Create empty electrode array
        electrodes = []

        # first coloumn
        elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14, rad_neck)  # x-coordinate of electrode
        elec_y = corr_y(dist_C7_Sp14, rad_neck)             # y-coordinate of electrode
        elec_z = loc_C7                                     # z-coordinate of electrode
        elec_pnt_OG = array([elec_x, elec_y, elec_z])       # coordinate of electrode in np.array-form for easier manipulation

        for i in range(5):
            elec_pnt = elec_pnt_OG + array([0,0,ver_dist_elec]) * i
            electrodes.append(elec_pnt)

        # second coloumn
        elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)  # x-coordinate of electrode
        elec_y = corr_y(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)             # y-coordinate of electrode
        elec_z = loc_C7 - rad_electrode                                             # z-coordinate of electrode
        elec_pnt_OG = array([elec_x, elec_y, elec_z])                               # coordinate of electrode in np.array-form for easier manipulation

        for i in range(6):
            elec_pnt = elec_pnt_OG + array([0,0,ver_dist_elec]) * i
            electrodes.append(elec_pnt)

        # third coloumn
        elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)  # x-coordinate of electrode
        elec_y = corr_y(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)             # y-coordinate of electrode
        elec_z = loc_C7                                                               # z-coordinate of electrode
        elec_pnt_OG = array([elec_x, elec_y, elec_z])                                 # coordinate of electrode in np.array-form for easier manipulation

        for i in range(5):
            elec_pnt = elec_pnt_OG + array([0,0,ver_dist_elec]) * i
            electrodes.append(elec_pnt)

        if show == True:
            import matplotlib.pyplot as plt
            
            corr_x = lambda d, r: r - r * cos( d / r )
            corr_y = lambda d,r: r * sin( d / r)

            rad_neck = params_geom[3][0]
            loc_neck = params_geom[3][1]
            dist_C7_Sp14 = 0.01
            rad_electrode = 0.009                   # radius electrode
            #ver_dist_elec = 0.018 



            circle = plt.Circle((loc_neck, 0), rad_neck, color = 'r', fill = False)
            _, ax = plt.subplots()

            ax.add_patch(circle)

            elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14, rad_neck)  # x-coordinate of electrode
            elec_y = corr_y(dist_C7_Sp14, rad_neck)             # y-coordinate of electrode

            elec_x2 = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)  # x-coordinate of electrode
            elec_y2 = corr_y(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)             # y-coordinate of electrode

            elec_x3 = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)  # x-coordinate of electrode
            elec_y3 = corr_y(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)             # y-coordinate of electrode

            ax.scatter(x = elec_x, y = elec_y)
            ax.scatter(x = elec_x2, y = elec_y2)
            ax.scatter(x = elec_x3, y = elec_y3)

            plt.xlim([-0.035, 0.115])
            plt.ylim([-0.075,0.075])
            #plt.show()

        print('Electrode positions calculated!')
        N_E = 16
        return electrodes, N_E

    elif electrode_pattern == '32':
        print('Alternatives are not yet available.')
        exit('Program closes now.')




