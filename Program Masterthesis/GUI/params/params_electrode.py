""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Material Conductivity values
    Author: Markus E. Oberndorfer
####################################################################### """


dist_C7_Sp14 = 0.01                     # horiz. dist. between C7 & Sp14 --> (Notation according to (Wimmer, Kostoglou, Müller-Putz; 2022) )
rad_electrode = 0.009                   # radius electrode
ver_dist_elec = 0.018                   # meter, vertical distance between electrodes

params_elec = [dist_C7_Sp14, rad_electrode, ver_dist_elec]