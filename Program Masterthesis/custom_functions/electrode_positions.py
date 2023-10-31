""" #######################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Calcualtion of Electrode Positions
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from math import sqrt, sin, cos
from numpy import array

# ----- FUNCTIONS -----
def ElectrodePositions(loc_C7, params_elec, params_geom, show = False):
    """
    Calculate the electrode positions for a given geometry and C7-location.
    The electrode pattern is as follwos:
        Three columns on the right side of the spine:
        One 16 elec. array consists of:
            Left Col.: 5 electrodes
            Middle Col.: 6 electrodes
            Right Col.: 5 electrodes
        Vertical distance between electrodes (center to center) is 18mm

    Input:
        - loc_C7: Location of the seventh cervical vertebral body. [float]
        - params_elec: Is a list of three floats that define the electrode
        setup. [list; float]
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
        - show: If set to True, a visualization of the electrod setup is
        shown. [bool]

    Output:
        - electrodes: is a list of arrays that contain the X,Y,Z coordinates
        of each electrode. [list, Numpy array]
        - N_E: is the total number of electrodes. [int]
    """
    dist_C7_Sp14 = params_elec[0]
    rad_electrode = params_elec[1]
    ver_dist_elec = params_elec[2]
    rad_neck = params_geom[3][0]
    loc_neck = params_geom[3][1]

    # Correction for curved neck surface
    corr_x = lambda d, r: r - r * cos( d / r )
    corr_y = lambda d,r: r * sin( d / r)

    # Create empty electrode array
    electrodes = []

    # first coloumn
    # x-coordinate of electrode
    elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14, rad_neck) 
    # y-coordinate of electrode
    elec_y = corr_y(dist_C7_Sp14, rad_neck)    
    # z-coordinate of electrode         
    elec_z = loc_C7                     
    # coordinate of electrode in np.array-form for easier manipulation                
    elec_pnt_OG = array([elec_x, elec_y, elec_z])       

    for i in range(5):
        elec_pnt = elec_pnt_OG + array([0,0,ver_dist_elec]) * i
        electrodes.append(elec_pnt)

    # second coloumn
    # x-coordinate of electrode
    elec_x = -(rad_neck - loc_neck) \
             + corr_x(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)
    # y-coordinate of electrode
    elec_y = corr_y(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)     
    # z-coordinate of electrode        
    elec_z = loc_C7 - rad_electrode                
    # coordinate of electrode in np.array-form for easier manipulation                             
    elec_pnt_OG = array([elec_x, elec_y, elec_z])                              

    for i in range(6):
        elec_pnt = elec_pnt_OG + array([0,0,ver_dist_elec]) * i
        electrodes.append(elec_pnt)

    # third coloumn
    # x-coordinate of electrode
    elec_x = -(rad_neck - loc_neck) \
             + corr_x(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)  
    # y-coordinate of electrode
    elec_y = corr_y(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)        
    # z-coordinate of electrode     
    elec_z = loc_C7  
    # coordinate of electrode in np.array-form for easier manipulation                                                             
    elec_pnt_OG = array([elec_x, elec_y, elec_z])                                 

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
        rad_electrode = 0.009

        circle = plt.Circle((loc_neck, 0), rad_neck, color = 'r',
                            fill = False)
        _, ax = plt.subplots()

        ax.add_patch(circle)
        # x-coordinate of electrode
        elec_x = -(rad_neck - loc_neck) + corr_x(dist_C7_Sp14, rad_neck)  
        # y-coordinate of electrode
        elec_y = corr_y(dist_C7_Sp14, rad_neck)             

        # x-coordinate of electrode
        elec_x2 = -(rad_neck - loc_neck) \
                  + corr_x(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)  
        # y-coordinate of electrode
        elec_y2 = corr_y(dist_C7_Sp14 + rad_electrode*sqrt(3), rad_neck)             

        # x-coordinate of electrode
        elec_x3 = -(rad_neck - loc_neck) \
                  + corr_x(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)
        # y-coordinate of electrode
        elec_y3 = corr_y(dist_C7_Sp14 + 2*rad_electrode*sqrt(3), rad_neck)             

        ax.scatter(x = elec_x, y = elec_y)
        ax.scatter(x = elec_x2, y = elec_y2)
        ax.scatter(x = elec_x3, y = elec_y3)

        plt.xlim([-0.035, 0.115])
        plt.ylim([-0.075,0.075])
        #plt.show()

    print('Electrode positions calculated!')
    N_E = 16
    return electrodes, N_E




