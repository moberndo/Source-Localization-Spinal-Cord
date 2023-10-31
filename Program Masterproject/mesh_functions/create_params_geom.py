from math import pi

def CreateParamsGeom(geom_str):
    geom_param = geom_str.rsplit(';')
    perimeter_neck = float(geom_param[0]) / 100
    Inion_C7_dist = float(geom_param[1]) / 100
    C7_L5_dist = float(geom_param[2]) / 100

    avg_neck_rad = 0.06

    

    ##### NECK ######
    rad_neck = perimeter_neck / (2*pi)      # 0.06
    loc_neck = rad_neck * 0.5
    len_neck = Inion_C7_dist * 0.9

    params_neck = [rad_neck, loc_neck, len_neck]

    ##### SPINAL CORD #####
    len_spine = len_neck + C7_L5_dist             # 0.7
    rad_spine = 0.01

    params_spine = [len_spine, rad_spine]

    ##### DISCS ######
    height_cervic_disc = 0.003
    height_thorac_disc = 0.005
    height_lumbar_disc = 0.006

    params_disc = [height_cervic_disc, height_thorac_disc, height_lumbar_disc]

    #### Cervical Vertebraes ####

    vert_cervic_height = (len_neck - 6 * height_cervic_disc) / 7
    # Atlas
    atlas_height = vert_cervic_height
    atlas_rad1 = ...
    atlas_rad2 = ...

    params_atlas = [atlas_height, atlas_rad1, atlas_rad2]

    # Axis
    axis_height = ...
    ...
    params_axis = [axis_height]

    # C3-C7
    C3C7_height = ...
    C3C7_rad_foramen = ...
    C3C7_rad_vert_body = ...

    params_C3C7 = [C3C7_height]
    # Combined

    params_vertebrae = [params_atlas, params_axis, params_C3C7]

    ##### COMBINED #####

    params_geom = [params_spine, params_vertebrae ,params_disc ,params_neck]#, params_trachea, params_oesoph]
    

    return params_geom

    '''
    

    

    ##### VERTEBRAE #####
    # assuming the outer radius of the vertical body is greater by a factor of 1.3 than then outer radius of the foramen, then:
    fac = 1.3
    rad_foramen = loc_neck / (0.85 * (fac + 1))              # 0.0135
    rad_vert_body = rad_foramen * fac                        # 0.0175

    # cervical part
    num_cervic_vert = 7
    height_cervic_vert = (Inion_C7_dist * 0.9 - height_cervic_disc * (num_cervic_vert  - 1)) / num_cervic_vert  # 0.015
    # thoracic part
    num_thorac_vert = 12
    height_thorac_vert =  (C7_L5_dist * 0.65 - height_thorac_disc * num_thorac_vert ) / num_thorac_vert            #   0.022
    # lumbar part
    num_lumbar_vert = 5
    height_lumbar_vert = (C7_L5_dist * 0.35 - height_lumbar_disc * num_lumbar_vert ) / num_lumbar_vert             # 0.025

    

    

    params_vertebrae = [rad_vert_body, rad_foramen, num_cervic_vert, height_cervic_vert, num_thorac_vert, 
                        height_thorac_vert, num_lumbar_vert, height_lumbar_vert]
    

    rad_oesoph = 0.008 *  rad_neck / avg_neck_rad
    rad_trachea = 0.009 *  rad_neck / avg_neck_rad

    dist = ( rad_neck - rad_vert_body - 2 * rad_trachea - 2 * rad_oesoph) / 3
    ##### TRACHEA ######
    
    loc_trachea = loc_neck + rad_vert_body + dist * 2 + rad_oesoph * 2 + rad_trachea                     #   0.075
    len_trachea = len_neck * 1.2            #  0.15
    thick_trachea = 0.002  *  rad_neck / avg_neck_rad

    params_trachea = [rad_trachea, loc_trachea, len_trachea, thick_trachea]

    ###### OESOPHAGUS ######
    
    loc_oesoph = loc_neck + rad_vert_body + dist + rad_oesoph               #   0.055
    len_oesoph = len_neck * 1.8             # 0.2
    thick_oesoph = 0.002 *  rad_neck / avg_neck_rad

    params_oesoph = [rad_oesoph, loc_oesoph, len_oesoph, thick_oesoph]
    '''