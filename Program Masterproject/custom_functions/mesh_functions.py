""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Functions to create mesh
    Author: Markus E. Oberndorfer
####################################################################### """

from ngsolve import Mesh, Draw
from netgen.csg import Cylinder, OrthoBrick, CSGeometry, Pnt
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
    len_neck = Inion_C7_dist

    params_neck = [rad_neck, loc_neck, len_neck]

    ##### DISCS ######
    height_cervic_disc = 0.003
    height_thorac_disc = 0.005
    height_lumbar_disc = 0.006

    params_disc = [height_cervic_disc, height_thorac_disc, height_lumbar_disc]

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

    ##### SPINAL CORD #####
    len_spine = (Inion_C7_dist + C7_L5_dist ) + 0.025            # 0.7
    rad_spine = 0.01

    params_spine = [len_spine, rad_spine]

    

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

    ##### COMBINED #####

    params_geom = [params_spine, params_vertebrae ,params_disc ,params_neck, params_trachea, params_oesoph]

    return params_geom

def create_mesh_1(length_spine, height_vertebral_body, radius_spine, rad_vert_body, rad_foramen, maxh_):
    height_disc = height_vertebral_body * 0.25
    height_unit = height_vertebral_body + height_disc
    start_pnt = 0

    def create_vertebre(height_vertebral_body, rad_foramen, start_pnt):
        vert_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), rad_foramen)
        vert = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_vertebral_body)) * vert_
        return vert

    def create_disc(rad_vert_body, start_pnt):  
        disc_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), rad_vert_body)
        disc = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_disc)) * disc_
        return disc

    def create_spine(length_spine, radius_spine):
        box1 = OrthoBrick(Pnt(-0.05,-0.05,-0.01), Pnt(0.05,0.05,length_spine-0.01))
        spine = Cylinder(Pnt(0,0,0), Pnt(0,0,length_spine), radius_spine) * box1
        return spine

    spine = create_spine(length_spine, radius_spine)
    vertebre = create_vertebre(height_vertebral_body, rad_foramen, start_pnt)
    disc = create_disc(rad_vert_body, start_pnt + height_vertebral_body)

    for i in range(1,33+1):
        start_point = height_unit * i
        vert = create_vertebre(height_vertebral_body, rad_foramen, start_point)
        vertebre = vertebre + vert
        start_point = height_vertebral_body + height_unit * i
        disc = disc + create_disc(rad_vert_body, start_point)

    vertebre = vertebre - spine
    disc = disc - spine
    
    geo = CSGeometry()
    geo.Add (spine.mat('spine')) 
    geo.Add (vertebre.mat('vertebre'))
    geo.Add (disc.mat('disc'))
    
    mesh = geo.GenerateMesh(maxh = maxh_)
    mesh.Save('./mesh_files/Spine_1.vol')
    mesh = Mesh(mesh)
    return mesh

def create_mesh_2(length_spine, height_vertebral_body, radius_spine, rad_vert_body, rad_foramen, maxh_):
    height_disc = height_vertebral_body * 0.25
    height_unit = height_vertebral_body + height_disc
    dist_vert_body = (rad_foramen + rad_vert_body)*0.85
    start_pnt = 0

    def create_vertebre(spine, height_vertebral_body, rad_vert_body, rad_foramen, start_pnt):
        #vertebre_1
        vert1_ = Cylinder(Pnt(dist_vert_body,0,0), Pnt(dist_vert_body,0,1), rad_vert_body)
        vert1 = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_vertebral_body)) * vert1_
        # vertebre_2
        vert2_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), rad_foramen)
        vert2 = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_vertebral_body)) * vert2_
        # vertebre_2_update
        remove_tube_ = OrthoBrick(Pnt(-1,-0.1,-0.1), Pnt(1,0.1,1))
        remove_tube1 = Cylinder(Pnt(-0.001,-1,start_pnt + height_vertebral_body+0.003), Pnt(-0.001,1,start_pnt + height_vertebral_body+0.003), 0.01) * remove_tube_
        remove_tube2 = Cylinder(Pnt(-0.001,-1,start_pnt-0.003), Pnt(-0.001,1,start_pnt-0.003), 0.01) * remove_tube_
        vert = vert1 + vert2
        vert = vert - spine
        vert = vert - remove_tube1
        vert = vert - remove_tube2
        return vert, remove_tube1, remove_tube2

    def create_disc(rad_vert_body, start_pnt, spine, remove_tube_1, remove_tube_2):  
        disc_ = Cylinder(Pnt(dist_vert_body,0,0), Pnt(dist_vert_body,0,1), rad_vert_body)
        disc = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_disc)) * disc_
        disc = disc - spine - remove_tube_1 - remove_tube_2
        return disc

    def create_spine(length_spine, radius_spine):
        box1 = OrthoBrick(Pnt(-0.05,-0.05,-0.01), Pnt(0.05,0.05,length_spine-0.01))
        spine = Cylinder(Pnt(0,0,0), Pnt(0,0,length_spine), radius_spine) * box1
        return spine

    spine = create_spine(length_spine, radius_spine)
    vertebre, remove_tube_1, remove_tube_2 = create_vertebre(spine, height_vertebral_body, rad_vert_body, rad_foramen, start_pnt)
    disc = create_disc(rad_vert_body, start_pnt + height_vertebral_body, spine, remove_tube_1, remove_tube_2)

    for i in range(1,33+1):
        start_point = height_unit * i
        vert, remove_tube_1, remove_tube_2 = create_vertebre(spine, height_vertebral_body, rad_vert_body, rad_foramen, start_point)
        vertebre = vertebre + vert
        start_point = height_vertebral_body + height_unit * i
        disc = disc + create_disc(rad_vert_body, start_point, spine, remove_tube_1, remove_tube_2)

    geo = CSGeometry()
    geo.Add (spine.mat('spine')) 
    geo.Add (vertebre.mat('vertebre'))
    geo.Add (disc.mat('disc'))
    
    mesh = geo.GenerateMesh(maxh = maxh_)
    mesh.Save('./mesh_files/Spine_2.vol')
    mesh = Mesh(mesh)
    return mesh

def create_mesh_3(params_spine, params_vert, params_disc, params_neck, params_oesoph, params_trachea,
                  geom_str, maxh_ = 0.005):

    def create_spine(length_spine, radius_spine):
        box1 = OrthoBrick(Pnt(-0.05,-0.05,-0.025), Pnt(0.05,0.05,length_spine + 0.025))
        spine = Cylinder(Pnt(0,0,0), Pnt(0,0,1), radius_spine) * box1
        return spine
    
    def create_vertebre(spine, start_pnt, rad_vert_body, rad_foramen, height_vertebral_body, dist_vert_body):
        #vertebre_1
        vert1_ = Cylinder(Pnt(dist_vert_body,0,0), Pnt(dist_vert_body,0,1), rad_vert_body)
        vert1 = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_vertebral_body)) * vert1_
        # vertebre_2
        vert2_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), rad_foramen)
        vert2 = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_vertebral_body)) * vert2_
        # vertebre_2_update
        remove_tube_ = OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1))
        rad_tube = height_vertebral_body / 3
        remove_tube1 = Cylinder(Pnt(0,-1,start_pnt + height_vertebral_body),
                                Pnt(0,1,start_pnt + height_vertebral_body), rad_tube * 0.5) * remove_tube_
        remove_tube2 = Cylinder(Pnt(0,-1,start_pnt), Pnt(0,1,start_pnt), rad_tube) * remove_tube_
        vert = vert1 + vert2
        vert = vert - spine
        vert = vert - remove_tube1
        vert = vert - remove_tube2
        return vert, remove_tube1, remove_tube2

    def create_disc(rad_vert_body, height_disc, start_pnt, spine, remove_tube_1, remove_tube_2, dist_vert_body): 
        disc_ = Cylinder(Pnt(dist_vert_body,0,0), Pnt(dist_vert_body,0,1), rad_vert_body)
        disc = OrthoBrick(Pnt(-1,-1,start_pnt), Pnt(1,1,start_pnt + height_disc)) * disc_
        disc = disc - spine - remove_tube_1 - remove_tube_2
        return disc
    
    def create_neck(params_neck, params_oesoph, params_trachea, start_neck):
    
        def create_trachea(rad_trachea, loc_trachea, len_trachea, thick_trachea, rad_neck, loc_neck, start_neck):
            start_trachea = Pnt(loc_trachea, 0, start_neck)
            end_trachea = Pnt(loc_trachea, 0, start_neck - len_trachea) 
            box_1 = OrthoBrick(Pnt(-rad_neck+loc_neck,-rad_neck,start_neck - len_trachea), Pnt(rad_neck+loc_neck,rad_neck,start_neck)) 
            trachea_outer = Cylinder(start_trachea, end_trachea, rad_trachea) * box_1
            trachea_inner = Cylinder(start_trachea, end_trachea, rad_trachea - thick_trachea) * box_1
            trachea = trachea_outer - trachea_inner
            return trachea, trachea_outer

        def create_oesophagus(rad_oesoph, loc_oesoph, len_oesoph, thick_oesoph, rad_neck, loc_neck, start_neck):
            start_oesophagus = Pnt(loc_oesoph, 0, start_neck)
            end_oesophagus = Pnt(loc_oesoph, 0, start_neck - len_oesoph)
            box_2 = OrthoBrick(Pnt(-rad_neck+loc_neck,-rad_neck,start_neck - len_oesoph), Pnt(rad_neck+loc_neck,rad_neck,start_neck)) 
            oesophagus_outer = Cylinder(start_oesophagus, end_oesophagus, rad_oesoph) * box_2
            oesophagus_inner = Cylinder(start_oesophagus, end_oesophagus, rad_oesoph-thick_oesoph) * box_2
            oesophagus = oesophagus_outer - oesophagus_inner
            return oesophagus, oesophagus_outer
        
        def create_neck_(rad_neck, loc_neck, len_neck, trachea_outer, oesoph_outer, start_neck):
            box_neck = OrthoBrick(Pnt(-rad_neck+loc_neck,-rad_neck,start_neck - len_neck), Pnt(rad_neck+loc_neck,rad_neck,start_neck)) 
            cyl_neck = Cylinder(Pnt(loc_neck, 0, -1), Pnt(loc_neck, 0, 1), rad_neck)
            neck = box_neck * cyl_neck - trachea_outer - oesoph_outer
            return neck
        
        trachea, trachea_outer = create_trachea(params_trachea[0], params_trachea[1], params_trachea[2], params_trachea[3],
                                                params_neck[0], params_neck[1], start_neck)
        oesophagus, oesophagus_outer = create_oesophagus(params_oesoph[0], params_oesoph[1], params_oesoph[2], params_oesoph[3],
                                                         params_neck[0], params_neck[1], start_neck)
        neck = create_neck_(params_neck[0], params_neck[1], params_neck[2], trachea_outer, oesophagus_outer, start_neck)
        
        return neck, trachea, oesophagus
    
    #spine_pnts = [] # are points where nerves exit the spinal cord

    spine = create_spine(params_spine[0], params_spine[1])
    dist_vert_body = ( params_vert[0] + params_vert[1] ) * 0.85

    # create cervical vertebraes    
    height_unit = params_vert[3] + params_disc[0]  
    start_vert_OG = params_spine[0] - height_unit         # z-coordinate of start point  
    for i in range(params_vert[2]):
        start_disc = start_vert_OG - height_unit * i
        start_vert = start_disc + params_disc[0]

        # create new vertebrae
        
        vert, remove_tube_1, remove_tube_2 = create_vertebre(spine, start_vert, params_vert[0], params_vert[1], params_vert[3], dist_vert_body)
        # check if vertebre already exists, if yes add the new one, if no, create a new variable/geometry
        if 'vertebre' in locals(): 
            vertebre = vertebre + vert
        else:
            vertebre = vert

        # create disc
        #start_disc = start_vert + params_vert[3]
        disc_ = create_disc(params_vert[0], params_disc[0], start_disc, spine, remove_tube_1, remove_tube_2, dist_vert_body)
        if 'disc' in locals():
            disc = disc + disc_
        else:
            disc = disc_


        
        
        
    loc_C7 = start_vert # - params_disc[0]

    # create thoracic vertebraes   
    height_unit = params_vert[5] + params_disc[1]  
    start_vert_OG = start_disc - height_unit         # z-coordinate of start point  
    for i in range(params_vert[4]):
        start_disc = start_vert_OG - height_unit * i
        start_vert = start_disc + params_disc[1]


        # create new vertebrae
        
        vert, remove_tube_1, remove_tube_2 = create_vertebre(spine, start_vert, params_vert[0], params_vert[1], params_vert[5], dist_vert_body)
        vertebre = vertebre + vert

        # create disc
        #start_vert = start_vert + params_vert[5]
        disc_ = create_disc(params_vert[0], params_disc[1], start_disc, spine, remove_tube_1, remove_tube_2, dist_vert_body)
        disc = disc + disc_

        
        

    # create lumbar vertebraes
    height_unit = params_vert[7] + params_disc[2]  
    start_vert_OG = start_disc - height_unit         # z-coordinate of start point  
    for i in range(params_vert[6]):
        start_disc = start_vert_OG - height_unit * i
        start_vert = start_disc + params_disc[2]


        # create new vertebrae
        
        vert, remove_tube_1, remove_tube_2 = create_vertebre(spine, start_vert, params_vert[0], params_vert[1], params_vert[7], dist_vert_body)
        vertebre = vertebre + vert

        # create disc
        #start_disc = start_vert + params_vert[7]
        disc_ = create_disc(params_vert[0], params_disc[2], start_disc, spine, remove_tube_1, remove_tube_2, dist_vert_body)
        disc = disc + disc_

        
        


    neck, trachea, oesophagus = create_neck(params_neck, params_oesoph, params_trachea, params_spine[0])
    neck = neck - (spine + disc + vertebre)

    geo = CSGeometry()
    geo.Add (spine.mat('spine')) 
    geo.Add (vertebre.mat('vertebre'))
    geo.Add (disc.mat('disc'))
    geo.Add (neck.mat('neck'))
    geo.Add (trachea.mat('trachea'))
    geo.Add (oesophagus.mat('oesophagus'))

    mesh = geo.GenerateMesh(maxh = maxh_)
    mesh.Save(geom_str[2])
    mesh = Mesh(mesh)
    # add new measurments_geom.txt file
    with open(geom_str[0], 'w') as f:
        f.write(geom_str[1])

    with open(geom_str[3], 'w') as f:
        f.write(str(loc_C7) + ';' + str(maxh_))

    return mesh.Curve(3), loc_C7, maxh_
    
if __name__ == '__main__':
    geom_str = '39;15;50'
    params_geom = CreateParamsGeom(geom_str)
    mesh = create_mesh_3(params_geom[0], params_geom[1], params_geom[2], params_geom[3],
                        params_geom[4], params_geom[5], geom_str)
    Draw(mesh)
