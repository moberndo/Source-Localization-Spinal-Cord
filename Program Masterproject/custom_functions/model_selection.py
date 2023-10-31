from pathlib import Path
from ngsolve import Mesh, Draw
from .mesh_functions import CreateParamsGeom, create_mesh_1, create_mesh_2, create_mesh_3

def ModelSelection(num_participant, main_path, show):

    # create functions for measurment check and read_data
    folder_participant = 'participant_' + num_participant
    #patha = Path(__file__)
    path_folder_participant = main_path / 'data' / folder_participant
    path_participant_geom = path_folder_participant / 'measurment_geom.txt'
    
    if path_participant_geom.exists():
        geom = open(path_participant_geom, 'r')
        measurments_geom = geom.read()
        geom.close()
    else:
        alt_path = main_path / 'data' / 'measurment_geom_alt.txt'
        if alt_path.exists():
            geom = open(alt_path, 'r')
            measurments_geom = geom.read()
            geom.close()
        else:
            exit('File missing! Alternative geometry measure not available.')

    params_geom = CreateParamsGeom(measurments_geom)


    path_vol_file = path_folder_participant / 'mesh_file.vol'
    path_spine_params = path_folder_participant / 'spine_params.txt'

    if path_vol_file.exists():
        mesh = Mesh(str(path_vol_file))
        mesh.Curve(3)
        if show == True:
            Draw(mesh)

        spine_params_ = open(path_spine_params, 'r')
        spine_params = spine_params_.read()
        spine_params_.close()

        spine_params = spine_params.rsplit(';')
        loc_C7 = float(spine_params[0])
        maxh_ = float(spine_params[1])

        print('All mesh files found and loaded.')

    else:
        print('Mesh of the newest model does not exist. Creating new Mesh now ...')
        mesh, loc_C7, maxh_ = create_mesh_3(params_geom[0], params_geom[1], params_geom[2], params_geom[3],
                             params_geom[4], params_geom[5],
                             [path_participant_geom, measurments_geom, path_vol_file, path_spine_params])
        print('Mesh sucessfully generated!')

    return mesh, params_geom, loc_C7, maxh_


if __name__ == '__main__':
    print('Everything loaded.')