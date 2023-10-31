""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to
    simulate the voltage distribution caused by a dipole
    Subtitle: Functions to create mesh
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from ngsolve import Mesh, Draw
from netgen.csg import Cylinder, OrthoBrick, CSGeometry, Pnt, Plane, Vec
from math import pi, atan

# ----- FUNCTIONS -----
def CreateParamsGeom(geom_str):
    """
    For a given geometry setting, this function calculates all
    necessary geoemtry parameters.

    Input:
        geom_str: String with three entries. See ReadMe.md file
        for more information. [str]

    Output: 
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
    """
    geom_param = geom_str.rsplit(';')
    perimeter_neck = float(geom_param[0]) / 100
    Inion_C7_dist = float(geom_param[1]) / 100
    C7_L5_dist = float(geom_param[2]) / 100

    avg_neck_rad = 0.06
    ratio_thoracic_lumbar = 0.6666

    ##### NECK ######
    neck_rad = perimeter_neck / (2*pi)    
    neck_dist = neck_rad * 0.5
    neck_height = Inion_C7_dist

    ##### SPINAL CORD #####
    spine_height = (Inion_C7_dist + C7_L5_dist ) 
    spine_rad = 0.01

    ##### DISCS ######
    # cervical
    cerv_disc_height = 0.003
    # throacic
    thoraric_disc_height = 0.005
    # lumbar
    lumbar_disc_height = 0.006

    ##### VERTEBRAE #####
    
    # cervical part
    cervic_vert_num = 7
    cervic_vert_rad_foramen_inner = spine_rad * 1.1
    cervic_vert_rad_foramen_outer = 0.0145
    cervic_vert_dist = cervic_vert_rad_foramen_outer * 0.95
    cervic_vert_rad_body = 0.008625  
    cervic_vert_height = (neck_height - 6 * cerv_disc_height) \
                        / cervic_vert_num
    cervic_vert_alpha_tubes = 20
    cervic_vert_rad_tubes = cervic_vert_rad_body * 1

    # axis
    axis_dist = cervic_vert_rad_foramen_outer * 0.95
    axis_rad_body = 0.008625
    axis_rad_foramen_inner = cervic_vert_rad_foramen_inner
    axis_rad_foramen_outer = cervic_vert_rad_foramen_outer
    axis_height = (neck_height - 6 * cerv_disc_height) / cervic_vert_num
    axis_alpha_tubes = 20
    axis_rad_tubes = axis_rad_body * 1
    axis_dens_rad = axis_rad_body * 0.85
    axis_dens_dist = axis_dist

    # atlas
    atlas_rad_inner = (axis_dist + axis_dens_rad \
                                 + axis_rad_foramen_inner) / 2
    atlas_dist = atlas_rad_inner - axis_rad_foramen_inner
    atlas_rad_outer = atlas_rad_inner * 1.25
    atlas_height = (neck_height - 6 * cerv_disc_height) / cervic_vert_num

    # thoracic part
    thorac_vert_num = 12
    thoracic_vert_rad_body = 0.012 
    thoracic_vert_dist = cervic_vert_dist - cervic_vert_rad_body \
                                          + thoracic_vert_rad_body
    
    thoracic_vert_rad_foramen_inner = spine_rad * 1.1
    thoracic_vert_rad_foramen_outer = cervic_vert_rad_foramen_outer
    thoracic_vert_height = (C7_L5_dist * ratio_thoracic_lumbar \
                            - thorac_vert_num * thoraric_disc_height) \
                            / thorac_vert_num  
    thoracic_vert_alpha_tubes = 15
    thoracic_vert_rad_tubes = thoracic_vert_rad_body * 1
    
    # lumbar part
    lumbar_vert_num = 5
    lumbar_vert_rad_body = 0.0125
    lumbar_vert_dist = cervic_vert_dist - cervic_vert_rad_body \
                     + lumbar_vert_rad_body
    lumbar_vert_rad_foramen_inner = spine_rad * 1.1 
    lumbar_vert_rad_foramen_outer = cervic_vert_rad_foramen_outer
    lumbar_vert_height = (C7_L5_dist * (1 - ratio_thoracic_lumbar) \
                       - lumbar_vert_num * lumbar_disc_height) \
                       / lumbar_vert_num
    lumbar_vert_alpha_tubes = 15
    lumbar_vert_rad_tubes = lumbar_vert_rad_body * 1

    ##### DISCS ######
    # cervical
    cerv_disc_dist = cervic_vert_rad_foramen_outer * 0.95
    cerv_dist_rad = cervic_vert_rad_body
    # throacic
    thoracic_disc_dist = thoracic_vert_dist
    thoracic_dist_rad = thoracic_vert_rad_body
    # lumbar
    lumbar_disc_dist = cervic_vert_dist - cervic_vert_rad_body \
                     + lumbar_vert_rad_body
    lumbar_disc_rad = lumbar_vert_rad_body

    ###### OESOPHAGUS ######
    oesoph_rad = 0.008 *  neck_rad / avg_neck_rad
    oesoph_height = neck_height * 1.8             
    oesoph_thick = 0.002 *  neck_rad / avg_neck_rad

    ##### TRACHEA ######
    trachea_rad = 0.009 *  neck_rad / avg_neck_rad
    trachea_height = neck_height * 1.2          
    trachea_thick = 0.002  *  neck_rad / avg_neck_rad

    ###### OESOPHAGUS & TRACHEA ######
    dist = ( neck_rad - cervic_vert_rad_body - 2 * trachea_rad - 2 \
            * oesoph_rad) / 3
    oesoph_dist = neck_dist + cervic_vert_rad_body + dist + oesoph_rad
    trachea_dist = neck_dist + cervic_vert_rad_body + dist * 2 \
                 + oesoph_rad * 2 + trachea_rad

    ##### COMBINED #####
    params_neck = [neck_rad, neck_dist, neck_height]

    params_spine = [spine_height, spine_rad]

    params_axis = [axis_dist, axis_rad_body, axis_rad_foramen_inner,
                   axis_rad_foramen_outer, axis_height, axis_alpha_tubes,
                   axis_rad_tubes,axis_dens_rad, axis_dens_dist]
    
    params_atlas = [atlas_dist, atlas_height, atlas_rad_inner,
                    atlas_rad_outer]

    params_cervic_vertebrae = [cervic_vert_num, cervic_vert_dist,
                               cervic_vert_rad_body,
                               cervic_vert_rad_foramen_inner,
                               cervic_vert_rad_foramen_outer,
                               cervic_vert_height,
                               cervic_vert_alpha_tubes, cervic_vert_rad_tubes]
    
    params_thoracic_vertebrae = [thorac_vert_num, thoracic_vert_dist,
                                 thoracic_vert_rad_body,
                                 thoracic_vert_rad_foramen_inner,
                                 thoracic_vert_rad_foramen_outer,
                                 thoracic_vert_height,
                                 thoracic_vert_alpha_tubes,
                                 thoracic_vert_rad_tubes]

    params_lumbar_vertebrae = [lumbar_vert_num, lumbar_vert_dist,
                               lumbar_vert_rad_body,
                               lumbar_vert_rad_foramen_inner,
                               lumbar_vert_rad_foramen_outer,
                               lumbar_vert_height,
                               lumbar_vert_alpha_tubes,
                               lumbar_vert_rad_tubes]

    params_vertebrae = [params_atlas, params_axis, params_cervic_vertebrae,
                        params_thoracic_vertebrae, params_lumbar_vertebrae]

    params_disc = [cerv_disc_dist, cerv_disc_height, cerv_dist_rad,
                   thoracic_disc_dist, thoraric_disc_height, 
                   thoracic_dist_rad, lumbar_disc_dist, lumbar_disc_height,
                   lumbar_disc_rad]
    
    params_trachea = [trachea_rad, trachea_dist, trachea_height,
                      trachea_thick]

    params_oesoph = [oesoph_rad, oesoph_dist, oesoph_height, oesoph_thick]
    
    params_geom = [params_spine, params_vertebrae ,params_disc ,params_neck,
                   params_trachea, params_oesoph]

    return params_geom

def CreateMesh(geom_str, save_file = None, maxh_ = 0.005):
    """
    This function creates the 3D model of the spine, including
    the spinal cord, the neck, the oesophagus and the trachea.

    Input:
        geom_str: String with three entries. See ReadMe.md file
        for more information. [str]
        - save_file: If this parameter is not none, it contains all
        path that are necessary to save the mesh.
        - maxh: Is the maximum size of one mesh element. [float]

    Output:
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - loc_C7: Location of the seventh cervical vertebral body. [float]
        - maxh_spine: Is the maximum size of one mesh element
        in the spinal cord. [float]
    """


    # All function-definitions that are following create the geometry that is
    # represented in their name. Therefore, they take some paramters of the 
    # parmas_geom variable and perform simple geometry manipulations like
    # the union or avergae quantity of two geoemtric bodies.
    def CreateNeck(start_pnt_neck, height_neck, rad_neck, dist_neck,
                   params_oesoph, params_trachea):
        def CreateOesophagus(start_pnt_oesoph, params_oesoph):
            oesoph_rad = params_oesoph[0]
            oesoph_dist = params_oesoph[1]
            oesoph_height = params_oesoph[2]
            oesoph_thick = params_oesoph[3]

            end_pnt_oesoph = start_pnt_oesoph - oesoph_height
            ortho_oesoph = OrthoBrick(Pnt(-100, -100, end_pnt_oesoph),
                                      Pnt(100, 100, start_pnt_oesoph))

            oesoph_outer = Cylinder(Pnt(oesoph_dist, 0, 0),
                                    Pnt(oesoph_dist, 0, 1),
                                    oesoph_rad) * ortho_oesoph
            oesoph_inner = Cylinder(Pnt(oesoph_dist, 0, 0),
                                    Pnt(oesoph_dist, 0, 1),
                                    oesoph_rad - oesoph_thick) * ortho_oesoph

            oesophagus = oesoph_outer - oesoph_inner

            return oesophagus, oesoph_outer

        def CreateTrachea(start_pnt_trachea, params_trachea):
            trachea_rad = params_trachea[0]
            trachea_dist = params_trachea[1]
            trachea_height = params_trachea[2]
            trachea_thick = params_trachea[3]

            end_pnt_trachea = start_pnt_trachea - trachea_height
            ortho_trachea = OrthoBrick(Pnt(-100, -100, end_pnt_trachea),
                                       Pnt(100, 100, start_pnt_trachea))

            trachea_outer = Cylinder(Pnt(trachea_dist, 0, 0),
                                     Pnt(trachea_dist, 0, 1),
                                     trachea_rad) * ortho_trachea
            trachea_inner = Cylinder(Pnt(trachea_dist, 0, 0),
                                     Pnt(trachea_dist, 0, 1),
                                     trachea_rad - trachea_thick) \
                                     * ortho_trachea

            trachea = trachea_outer - trachea_inner

            return trachea, trachea_outer

        end_pnt_neck = start_pnt_neck - height_neck  * 1.1
        ortho_neck = OrthoBrick(Pnt(-100, -100, end_pnt_neck),
                                Pnt(100, 100, start_pnt_neck))
        neck = Cylinder(Pnt(dist_neck, 0, 0),
                        Pnt(dist_neck, 0, 1), rad_neck) * ortho_neck
        trachea, trachea_outer = CreateTrachea(start_pnt_neck,
                                               params_trachea)
        oesophagus, oesophagus_outer = CreateOesophagus(start_pnt_neck,
                                                        params_oesoph)
        neck = neck - trachea_outer - oesophagus_outer
        return neck, trachea, oesophagus

    def CreateSpine(start_pnt_spine, height_spine, rad_spine):
        maxh_spine = 0.001
        end_pnt_spine = start_pnt_spine - height_spine
        ortho_spine = OrthoBrick(Pnt(-100, -100, end_pnt_spine),
                                 Pnt(100, 100, start_pnt_spine))
        spine = Cylinder(Pnt(0,0,0), Pnt(0,0,1),
                         rad_spine*1.01).maxh(maxh_spine) * ortho_spine
        return spine, maxh_spine

    def CreateAtlas(start_pnt_atlas, height_atlas, atlas_dist,
                    atlas_rad_outer, atlas_rad_inner):
        end_pnt_atlas = start_pnt_atlas - height_atlas
        ortho_atlas = OrthoBrick(Pnt(-100, -100, end_pnt_atlas \
                                     + height_atlas * 0.05),
                                 Pnt(100, 100, start_pnt_atlas))
        atlas_outer = Cylinder(Pnt(atlas_dist, 0, 0),
                               Pnt(atlas_dist, 0, 1),
                               atlas_rad_outer) * ortho_atlas
        atlas_inner = Cylinder(Pnt(atlas_dist, 0, 0),
                               Pnt(atlas_dist, 0, 1),
                               atlas_rad_inner) * ortho_atlas
        atlas = atlas_outer - atlas_inner
        return atlas, end_pnt_atlas

    def CreateAxis(end_pnt_atlas, height_axis, axis_dens_rad, axis_dist_spine,
                   rad_body,
                rad_foramen_outer, rad_foramen_inner, alpha_tubes, rad_tubes):
        axis, _ = CreateVertebrae(end_pnt_atlas,height_axis, axis_dist_spine,
                                  rad_body, rad_foramen_outer,
                                  rad_foramen_inner,
                            alpha_tubes, rad_tubes)
        start_pnt_axis = end_pnt_atlas
        end_pnt_axis = start_pnt_axis - height_axis
        # -- Update: Dens --
        axis_dens_rad_ = axis_dens_rad * 0.95
        dens_ = Cylinder(Pnt(axis_dist_spine, 0, 0),
                         Pnt(axis_dist_spine, 0, 1),
                         axis_dens_rad_)
        dens_brick = OrthoBrick(Pnt(-100, -100, end_pnt_axis),
                                Pnt(100, 100, start_pnt_axis + height_axis))
        dens = dens_ * dens_brick
        axis = axis + dens
        return axis, end_pnt_axis

    def CreateDisc(start_pnt_disc, height_disc, dist_disc, rad_disc, type, i):
        if i == 0 and (type == 'thoracic' or type == 'lumbar'):
            end_pnt_disc = start_pnt_disc - height_disc
            start_pnt_disc_ = start_pnt_disc * 0.99935
            ortho_disc = OrthoBrick(Pnt(-100, -100, end_pnt_disc),
                                    Pnt(100,100, start_pnt_disc_))
        else:
            end_pnt_disc = start_pnt_disc - height_disc
            ortho_disc = OrthoBrick(Pnt(-100, -100, end_pnt_disc),
                                    Pnt(100, 100, start_pnt_disc))
        disc = Cylinder(Pnt(dist_disc, 0, 0),
                        Pnt(dist_disc, 0, 1), rad_disc) * ortho_disc
        return disc, end_pnt_disc

    def CreateVertebrae(end_pnt_disc, height_vert, dist_body, rad_body,
                        rad_foramen_outer, rad_foramen_inner, alpha_tubes,
                        rad_tubes):
        start_pnt_vert = end_pnt_disc
        end_pnt_vert = start_pnt_vert - height_vert
        ortho_vert = OrthoBrick(Pnt(-100, -100, end_pnt_vert),
                                Pnt(100, 100, start_pnt_vert))

        # -- Body --
        body = Cylinder(Pnt(dist_body, 0, 0), Pnt(dist_body, 0, 1),
                        rad_body) * ortho_vert

        # -- Foramen --
        foramen_outer = Cylinder(Pnt(0,0,0), Pnt(0,0,1),
                                 rad_foramen_outer) * ortho_vert
        foramen_inner = Cylinder(Pnt(0,0,0), Pnt(0,0,1),
                                 rad_foramen_inner) * ortho_vert
        foramen = foramen_outer - foramen_inner
        vert = foramen + body

        # -- Update: Spinous Posterior --
        const1 = rad_foramen_outer / 2
        plane1 = Plane(Pnt(0,0,0), Vec(1,0,0))
        plane2 = Plane(Pnt(-const1*2, const1 , 0), Vec(-0.5,1,0))
        plane3 = Plane(Pnt(-const1*2,-const1, 0), Vec(-0.5,-1,0))
        cylinder1 = Cylinder(Pnt(-const1*3.25, -const1*2.25,0),
                             Pnt(-const1*3.25, -const1*2.25,1),
                             rad_foramen_outer*1.0) * ortho_vert
        cylinder2 = Cylinder(Pnt(-const1*3.25, const1*2.25,0),
                             Pnt(-const1*3.25, const1*2.25,1),
                             rad_foramen_outer*1.0) * ortho_vert
        spinous_post = plane1 * plane2 * plane3 * ortho_vert \
                     - cylinder1 - cylinder2 - foramen_inner
        plane4 = Plane(Pnt(0,0,start_pnt_vert),
                       Vec(-0.1, 0, 0.8)) \
                       * OrthoBrick(Pnt(-0.05, -0.05, end_pnt_vert*0.95),
                                    Pnt(0.05, 0.05, start_pnt_vert *1.05))
        plane5 = Plane(Pnt(0,0,end_pnt_vert + height_vert * 0.3),
                       Vec(0.1, 0, -0.8)) \
                       * OrthoBrick(Pnt(-0.05, -0.05, end_pnt_vert*0.95),
                                    Pnt(0.05, 0.05, start_pnt_vert *1.05))
        plane5 = plane5 * plane1
        spinous_post = spinous_post * plane4 * plane5
        vert = vert + spinous_post
        
        # -- Update: Remove Tubes --
        alpha = alpha_tubes * pi / 180 
        x_comp = atan(alpha)
        const2 = -rad_foramen_inner*0.2
        plane1 = Plane(Pnt(0,0,start_pnt_vert), Vec(-x_comp,-1, 0))
        plane2 = Plane(Pnt(0,0,start_pnt_vert), Vec(-x_comp,1, 0))
        tube1_upper = Cylinder(Pnt(const2,0,start_pnt_vert),
                               Pnt(x_comp, 1, start_pnt_vert),
                               rad_tubes*0.4) * plane1
        tube2_upper = Cylinder(Pnt(const2,0,start_pnt_vert), 
                               Pnt(x_comp, -1, start_pnt_vert),
                               rad_tubes*0.4) * plane2
        tube1_lower = Cylinder(Pnt(const2,0,end_pnt_vert),
                               Pnt(x_comp, 1, end_pnt_vert),
                               rad_tubes*0.4) * plane1
        tube2_lower = Cylinder(Pnt(const2,0,end_pnt_vert),
                               Pnt(x_comp, -1, end_pnt_vert),
                               rad_tubes*0.4) * plane2
        vert = vert - tube1_upper - tube1_lower - tube2_upper - tube2_lower

        return vert, end_pnt_vert

    params_geom = CreateParamsGeom(geom_str)
    # spine
    spine_height = params_geom[0][0] * 1
    spine_rad = params_geom[0][1]  * 1
    # atlas
    atlas_dist = params_geom[1][0][0] * 1
    atlas_rad_inner = params_geom[1][0][2] * 1
    atlas_rad_outer = params_geom[1][0][3] * 1
    atlas_height = params_geom[1][0][1] * 1
    # axis
    axis_dist = params_geom[1][1][0] * 1
    axis_rad_body = params_geom[1][1][1] * 1
    axis_rad_foramen_inner = params_geom[1][1][2] * 1
    axis_rad_foramen_outer = params_geom[1][1][3] * 1
    axis_height = params_geom[1][1][4] * 1
    axis_alpha_tubes = params_geom[1][1][5]
    axis_rad_tubes = params_geom[1][1][6] * 1
    axis_dens_rad = params_geom[1][1][7] * 1
    axis_dens_dist = params_geom[1][1][8] * 1
    # cervical disc
    cervic_disc_dist = params_geom[2][0] * 1
    cervic_disc_height = params_geom[2][1] * 1
    cervic_disc_rad = params_geom[2][2] * 1
    # cervical vertebrae
    cervic_vert_dist = params_geom[1][2][1] * 1
    cervic_vert_rad_body = params_geom[1][2][2] * 1
    cervic_vert_rad_foramen_inner = params_geom[1][2][3] * 1
    cervic_vert_rad_foramen_outer = params_geom[1][2][4] * 1
    cervic_vert_height = params_geom[1][2][5] * 1
    cervic_vert_alpha_tubes = params_geom[1][2][6]
    cervic_vert_rad_tubes = params_geom[1][2][7] * 1
    # thoracic disc
    thoracic_disc_dist = params_geom[2][3] * 1
    thoracic_disc_height = params_geom[2][4] * 1
    thoracic_disc_rad = params_geom[2][5] * 1
    # thoracic vertebrae
    thoracic_vert_dist = params_geom[1][3][1] * 1
    thoracic_vert_rad_body = params_geom[1][3][2] * 1
    thoracic_vert_rad_foramen_inner = params_geom[1][3][3] * 1
    thoracic_vert_rad_foramen_outer = params_geom[1][3][4] * 1
    thoracic_vert_height = params_geom[1][3][5] * 1
    thoracic_vert_alpha_tubes = params_geom[1][3][6]
    thoracic_vert_rad_tubes = params_geom[1][3][7] * 1
    # lumbar disc
    lumbar_disc_dist = params_geom[2][6] * 1
    lumbar_disc_height = params_geom[2][7] * 1
    lumbar_disc_rad = params_geom[2][8] * 1
    # lumbar vertebrae
    lumbar_vert_dist = params_geom[1][4][1] * 1
    lumbar_vert_rad_body = params_geom[1][4][2] * 1
    lumbar_vert_rad_foramen_inner = params_geom[1][4][3] * 1
    lumbar_vert_rad_foramen_outer = params_geom[1][4][4] * 1
    lumbar_vert_height = params_geom[1][4][5] * 1
    lumbar_vert_alpha_tubes = params_geom[1][4][6]
    lumbar_vert_rad_tubes = params_geom[1][4][7] * 1
    # neck
    neck_rad = params_geom[3][0] * 1
    neck_dist = params_geom[3][1] * 1
    neck_height = params_geom[3][2] * 1
    # oesophagus
    params_oesoph = params_geom[5]
    # trachea
    params_trachea = params_geom[4]
    
    neck, trachea, oesophaugs = CreateNeck(spine_height, neck_height,
                                           neck_rad, neck_dist,
                                           params_oesoph, params_trachea)
    spine, maxh_spine = CreateSpine(spine_height, spine_height, spine_rad)
    atlas, end_pnt_atlas = CreateAtlas(spine_height, atlas_height, atlas_dist,
                                       atlas_rad_outer, atlas_rad_inner)
    axis, start_pnt_og = CreateAxis(end_pnt_atlas, axis_height, axis_dens_rad,
                                    axis_dist, axis_rad_body,
                                    axis_rad_foramen_outer,
                                    axis_rad_foramen_inner, axis_alpha_tubes,
                                    axis_rad_tubes)
    vertebraes = atlas + axis 
    for i in range(5):
        start_pnt_disc = start_pnt_og - i \
                       * (cervic_disc_height + cervic_vert_height)
        disc_i, end_pnt_disc_i = CreateDisc(start_pnt_disc,
                                            cervic_disc_height,
                                            cervic_disc_dist, cervic_disc_rad,
                                            'cervical', i)
        out_ = CreateVertebrae(end_pnt_disc_i,cervic_vert_height,
                               cervic_vert_dist, cervic_vert_rad_body,
                               cervic_vert_rad_foramen_outer,
                               cervic_vert_rad_foramen_inner,
                               cervic_vert_alpha_tubes,
                                cervic_vert_rad_tubes)
        
        vertebrae_i = out_[0]
        end_pnt_vert_i = out_[1]
        
        vertebraes = vertebraes + vertebrae_i
        if 'discs' in locals(): 
            discs = discs + disc_i
        else:
            discs = disc_i

    loc_C7 = end_pnt_vert_i
    start_pnt_og = start_pnt_og - 5 \
                 * (cervic_disc_height + cervic_vert_height)

    for i in range(12):
        start_pnt_disc = start_pnt_og - i \
                       * (thoracic_disc_height + thoracic_vert_height)
        disc_i, end_pnt_disc_i = CreateDisc(start_pnt_disc,
                                            thoracic_disc_height,
                                            thoracic_disc_dist,
                                            thoracic_disc_rad,
                                            'thoracic', i)
        vertebrae_i, _ = CreateVertebrae(end_pnt_disc_i, thoracic_vert_height,
                                         thoracic_vert_dist,
                                         thoracic_vert_rad_body,
                                         thoracic_vert_rad_foramen_outer,
                                         thoracic_vert_rad_foramen_inner,
                                         thoracic_vert_alpha_tubes,
                                         thoracic_vert_rad_tubes)
        discs = discs + disc_i
        vertebraes = vertebraes + vertebrae_i

    
    start_pnt_og = start_pnt_og - 12 \
                 * (thoracic_disc_height + thoracic_vert_height)
    for i in range(5):
        start_pnt_disc = start_pnt_og - i \
                       * (lumbar_disc_height + lumbar_vert_height)
        disc_i, end_pnt_disc_i = CreateDisc(start_pnt_disc,
                                            lumbar_disc_height,
                                            lumbar_disc_dist,
                                            lumbar_disc_rad,
                                            'lumbar', i)
        vertebrae_i, _ = CreateVertebrae(end_pnt_disc_i, lumbar_vert_height,
                                         lumbar_vert_dist,
                                         lumbar_vert_rad_body,
                                         lumbar_vert_rad_foramen_outer, 
                                         lumbar_vert_rad_foramen_inner,
                                         lumbar_vert_alpha_tubes,
                                         lumbar_vert_rad_tubes)
        
        discs = discs + disc_i
        vertebraes = vertebraes + vertebrae_i
        
    discs = discs - spine
    vertebraes = vertebraes - spine
    neck = neck - spine - discs - vertebraes

    geo = CSGeometry()

    geo.Add(discs.mat('cartilage').col([1,0,0]))
    geo.Add(vertebraes.mat('bone').col([0,0,1]))
    geo.Add(spine.mat('nervous_tissue').col([0,1,0]))
    geo.Add(neck.mat('muscle').col([0.5,0,0.5]))
    geo.Add(oesophaugs.mat('oesophagus').col([0,0.2,0.8])) 
    geo.Add(trachea.mat('trachea').col([0.8,0,0.2]))
    
    geo.Draw()
    
    mesh = geo.GenerateMesh(maxh = maxh_)

    if save_file != None:
        mesh.Save(save_file[0])
        mesh = Mesh(mesh)
        # add new measurments_geom.txt file
        with open(save_file[1], 'w') as f:
            f.write(geom_str)

        with open(save_file[2], 'w') as f:
            f.write(str(loc_C7) + ';' + str(maxh_))
    else:
        mesh = Mesh(mesh)

    return mesh.Curve(3), loc_C7, maxh_spine


if __name__ == '__main__':
    # imports
    from ngsolve import Parameter
    from ngsolve import ElementId, VOL, TaskManager, GridFunction, CGSolver
    from ngsolve import H1, CoefficientFunction, BilinearForm, LinearForm, grad, dx, Preconditioner
    from math import sin, cos
    # constants
    maxh_ = 0.005
    geom_str = '39;15;50'
    mat_conduc = {"nervous_tissue" : 0.22, 
                "bone" : 0.0014, 
                "cartilage" : 0.008, 
                "muscle" : 0.01 , 
                "trachea" : 0.015 ,
                "oesophagus" : 0.015}
    # functions
    def FEM_Params(mesh, conductivity_vals):
        fes = H1(mesh, order = '1')     # define Finite Element Space
        u,v = fes.TnT()                 # define Trial and Test Space

        sigma_coeff = CoefficientFunction([conductivity_vals[mat] for mat in mesh.GetMaterials()])

        # create bounded X-elliptic bilinear form
        a = BilinearForm(fes)
        a += grad(v)*sigma_coeff*grad(u) * dx + 1e-8*sigma_coeff*u*v*dx
        c = Preconditioner(a, "bddc")

        # create bounded linear form
        f = LinearForm(fes)
        print('FEM ready!')

        return a, c, f, fes

    def MakeDipole3D(f, x, y, z, max_mesh, alpha, beta, p):
        spc = f.space
        mp1 = spc.mesh(x,y,z)
        ei1 = ElementId(VOL, mp1.nr)    
        fel1 = spc.GetFE(ei1)
        dnums1 = spc.GetDofNrs(ei1)
        shape1 = fel1.CalcShape(*mp1.pnt)

        for d1,s1 in zip(dnums1, shape1):
            f.vec[d1] += ( p.Get()*s1 )

        x1 = x + max_mesh * cos(beta) * cos(alpha)
        y1 = y + max_mesh * cos(beta) * sin(alpha)
        z1 = z + max_mesh * sin(beta)

        if x1 > 0.01:
            x1 = 0.01
            print('x1 too large, cap at 0.01.')
        if y1 > 0.01:
            print(y1)
            y1 = 0.01
            print('y1 too large, cap at 0.01.')

        mp2 = spc.mesh(x1, y1, z1)
        ei2 = ElementId(VOL, mp2.nr)
        fel2 = spc.GetFE(ei2)
        dnums2 = spc.GetDofNrs(ei2)
        shape2 = fel2.CalcShape(*mp2.pnt)

        for d2, s2 in zip(dnums2, shape2):
            f.vec[d2] += ( - p.Get()*s2 )

    def CalcGFU(dipoles,f,maxh_, a, c, fes):
        for dip in dipoles: 
            x = dip[0]
            y = dip[1]
            z = dip[2]
            alpha = dip[3]
            beta = dip[4]
            p = dip[5]
            MakeDipole3D(f, x = x, y = y, z = z, max_mesh= maxh_*1.005, alpha= alpha, beta=beta, p = p)

        # Assemble equations
        with TaskManager():
            a.Assemble()

        # Solve PDE
        gfu = GridFunction(fes)
        inv = CGSolver(a.mat, c.mat, printrates=False, precision=1e-8)
        gfu.vec.data = inv * f.vec

        return gfu, inv, f.vec
    # main
    print('succesfull start')
    params_geom = CreateParamsGeom(geom_str=geom_str)
    spine_len = params_geom[0][0]
    spine_rad = params_geom[0][1]
    neck_height = params_geom[3][2]

    print('starting with mesh now')
    mesh, _, _ = CreateMesh(geom_str, save_file=None)
    print('mesh created')
    Draw(mesh)
    #sigma_coeff = CoefficientFunction([mat_conduc[mat] for mat in mesh.GetMaterials()])
    #Draw(mesh)
    #Draw(sigma_coeff, mesh, 'sigma_coeff', draw_surf = True)