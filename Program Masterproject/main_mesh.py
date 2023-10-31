from ngsolve import Mesh, Draw
from math import pi
from netgen.csg import Cylinder, OrthoBrick, CSGeometry, Pnt, Plane, EllipticCone, Vec
from netgen.csg import CSGeometry

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

    # ==== Atlas ====
    # height
    atlas_height = vert_cervic_height
    # foramen
    atlas_rad_foramen_outer = 0.0175
    atlas_rad_foramen_inner = 0.0155
    # body
    atlas_rad1 = 0.0135
    atlas_rad2 = atlas_rad1 * 0.5
    atlas_dist_spine = atlas_rad_foramen_outer * 0.75
    # -- combined --
    params_atlas = [atlas_height, atlas_rad1, atlas_rad2]

    # ==== Axis ====
    # height
    axis_height = vert_cervic_height
    # foramen
    axis_rad_foramen_outer = 0.0175
    axis_rad_foramen_inner = 0.0145
    # body
    axis_rad1 = 0.0135
    axis_rad2 = axis_rad1 * 0.55
    axis_dist_spine = axis_rad_foramen_outer * 0.75
    # dens
    axis_dens_rad = axis_rad2 * 0.85
    axis_dens_height = axis_height * 0.8
    # -- combined --
    params_axis = [axis_height, axis_rad_foramen_outer, axis_rad_foramen_inner, axis_rad1, axis_rad2,
                   axis_dist_spine, axis_dens_rad, axis_dens_height]


    # ==== C3-C7 ====
    C3C7_height = ...
    C3C7_rad_foramen = ...
    C3C7_rad_vert_body = ...

    params_C3C7 = [C3C7_height]
    # Combined

    params_vertebrae = [params_atlas, params_axis, params_C3C7]

    ##### COMBINED #####

    params_geom = [params_spine, params_vertebrae ,params_disc ,params_neck]#, params_trachea, params_oesoph]
    

    return params_geom

def CreateAtlas(params_geom):
    # 
    len_spine = params_geom[0][0]
    height_atlas = params_geom[1][0][0]

    start_atlas = len_spine
    end_atlas = len_spine - height_atlas

    atlas_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), 0.2)
    atlas = OrthoBrick(Pnt(-1,-1,end_atlas), Pnt(1,1,start_atlas)) * atlas_
    return atlas, end_atlas

def CreateAxis(params_geom, start_axis):
    axis_height = params_geom[1][1][0]
    axis_rad_foramen_outer = params_geom[1][1][1]
    axis_rad_foramen_inner = params_geom[1][1][2]
    axis_rad1 = params_geom[1][1][3]
    axis_rad2 = params_geom[1][1][4]
    axis_dist_spine = params_geom[1][1][5]
    axis_dens_rad = params_geom[1][1][6]
    axis_dens_height = params_geom[1][1][7]

    
    end_axis = start_axis - axis_height
    axis_foramen_brick = OrthoBrick(Pnt(-1,-1,end_axis), Pnt(1,1,start_axis))

    # spinous posterior

    plane1 = Plane(Pnt(0,0,0), Vec())

    # foramen
    
    axis_foramen_inner = Cylinder(Pnt(0,0,0), Pnt(0,0,1), axis_rad_foramen_inner) * axis_foramen_brick
    axis_formane_outer = Cylinder(Pnt(0,0,0), Pnt(0,0,1), axis_rad_foramen_outer) * axis_foramen_brick
    
    axis_foramen =  axis_formane_outer - axis_foramen_inner

    # body
    axis_body = EllipticCone(a=Pnt(axis_dist_spine, 0, end_axis), vl=Vec(axis_rad2,0,0), vs=Vec(0,axis_rad1,0),
                             h= end_axis + axis_height, r=1) \
                * Plane(Pnt(0,0,end_axis), Vec(0, 0, -1)) * Plane(Pnt(0,0,start_axis), Vec(0,0,1))

    # dens
    dens_ = Cylinder(Pnt(axis_dist_spine, 0, 0), Pnt(axis_dist_spine, 0, 1), axis_dens_rad)
    dens_brick = OrthoBrick(Pnt(-1, -1, end_axis), Pnt(1, 1, start_axis + axis_dens_height*0.8))
    dens_top = EllipticCone(a=Pnt(axis_dist_spine, 0, start_axis + axis_dens_height*0.8), vl=Vec(axis_dens_rad,0,0),
                            vs=Vec(0,axis_dens_rad,0), h=axis_dens_height*0.2, r=0.5) \
                * Plane(Pnt(0,0,start_axis + axis_dens_height*0.8), Vec(0, 0, -1)) * Plane(Pnt(0,0,start_axis + axis_dens_height), Vec(0,0,1))
    dens = dens_ * dens_brick + dens_top

    
    



    # add geometries
    axis = axis_foramen + axis_body + dens

    return axis, end_axis

def CreateC3C7():
    ...
    return C3C7, end_C3C7

def CreateVertDiscs():
    ...


def CreateCervicalSpine(params_geom, geo):
    atlas, end_atlas = CreateAtlas(params_geom)
    axis, end_axis = CreateAxis(params_geom, end_atlas)

    cervical_spine = axis

    geo.Add(cervical_spine.mat('bone'))
    return geo


geom_str = '39;15;50'
params_geom = CreateParamsGeom(geom_str)

geo = CSGeometry()
geo = CreateCervicalSpine(params_geom, geo)

mesh = geo.GenerateMesh(maxh = 0.005)
mesh = Mesh(mesh)
Draw(mesh)

