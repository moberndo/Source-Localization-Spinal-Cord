"""####################################################################
    Title: Solving the source localization problem of spinal cord
    potentials for electric dipoles in the cervical spinal cord
    Subtitle: Select a 3D model
    Author: Markus E. Oberndorfer
####################################################################"""

# ----- IMPORT PACKAGES -----
from ngsolve import Mesh

# ----- IMPORT CUSTOM FUNCTIONS -----
from .mesh_functions import CreateParamsGeom, CreateMesh

# ----- FUNCTIONS -----
def ModelSelection(measurments_geom, main_path, mesh_folder):
    """
    This function either load an existing mesh, or created a new one.

    Input: 
        - user_ID: Is a short ID that points to the folder structure of
        interest. [str]
        - main_path: Contains the path of the main file [pathlib object]
    
    Output:
        - mesh: Is an object that stores the model geometry and its
        properties. [NGSolve mesh object]
        - params_geom: Software intern structure. See ReadMe.txt file
          for more information. [-]
        - loc_C7: Location of the seventh cervical vertebral body. [float]
        - maxh_: Is the maximum size of one mesh element. [float]
    """

    # create functions for measurment check and read_data
    params_geom = CreateParamsGeom(measurments_geom)

    path_vol_file = mesh_folder / 'mesh_file.vol'
    path_spine_params = mesh_folder / 'spine_params.txt'

    if path_vol_file.exists():
        mesh = Mesh(str(path_vol_file))
        mesh.Curve(3)

        spine_params_ = open(path_spine_params, 'r')
        spine_params = spine_params_.read()
        spine_params_.close()

        spine_params = spine_params.rsplit(';')
        loc_C7 = float(spine_params[0])
        maxh_ = float(spine_params[1])

        print('All mesh files found and loaded.')
    else:
        print('Mesh of the newest model does not exist. \
              Creating new Mesh now ...')
        save_file = [path_vol_file, path_spine_params]
        mesh, loc_C7, maxh_ = CreateMesh(measurments_geom, save_file)
        print('Mesh sucessfully generated!')

    return mesh, params_geom, loc_C7, maxh_


if __name__ == '__main__':
    ...