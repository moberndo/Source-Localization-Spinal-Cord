from ngsolve import Mesh, Draw
from netgen.csg import Cylinder, OrthoBrick, CSGeometry, Pnt
from math import pi

def CreateAtlas(params_geom):
    # 
    len_spine = params_geom[0][0]
    height_atlas = params_geom[1][0][0]

    start_atlas = len_spine
    end_atlas = len_spine - height_atlas

    atlas_ = Cylinder(Pnt(0,0,0), Pnt(0,0,1), 0.2)
    atlas = OrthoBrick(Pnt(-1,-1,start_atlas), Pnt(1,1,end_atlas)) * atlas_
    return atlas, end_atlas

def CreateAxis():
    ...
    return axis, end_axis

def CreateC3C7():
    ...
    return C3C7, end_C3C7

def CreateVertDiscs():
    ...


def CreateVertebraelSpine(params_geom, geo):
    atlas, end_atlas = CreateAtlas(params_geom)


    geo.Add(atlas)
    return geo
    

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
