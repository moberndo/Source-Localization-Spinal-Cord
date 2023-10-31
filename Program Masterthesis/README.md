# Solving the source localization problem of spinal cord potentials for electric dipoles in the cervical spinal cord

This is the software package related to the master thesis:
*__Solving the source localization problem of spinal cord potentials for electric dipoles in the cervical spinal cord - Markus E. Oberndorfer__*. 

**Version 1.0 \
November, 2023**

## Motivation

As the title states, the goal of this code is to finde the sources from where the spinal cord potentials originate.
This work is one of the first attempts to solve this kind of inverse problem.

## General project structure

The general structure of this project is as follows:

    - data
        - measurement_geom_alt.txt
        + user_XYZ
    - main.py
    + custom_functions
    + params

In the following, all folders with their respective contents are described.

### data - user_XYZ

Contains the spinal cord potential measurments per participant in the folder SpinalCordPotentials as well as the processed data.
The SpinalCordPotentials folder is organized per participant, in which each participant is expected to have several runs.
In the processed_data folder, the overall average is stored as data_avg.npy and the average per participant is stored as data_avg_part.npy.

    - data
        - measurement_geom_alt.txt
        - user_XYZ
            + images
            - measurement_geom.txt
            - SpinalCordPotentials
                - participant_1
                    - run_1.mat
                    - run_2.mat
                    - ...
                - participant_2
                    - run_1.mat
                    - run_2.mat
                    - ...
                - ...
            - processed_data
                - data_avg.npy
                - data_avg_part.npy
        - user_UVW
            + images
            - measurement_geom.txt
            - ...


### custom_functions

All custom functions needed for the main file are stored here. To get information on each function inside those files use the help(...) command.

    - custum_functions
        - dipole_positions.py
        - electrode_positions.py
        - evaluate.py
        - fem_functions.py
        - inv_functions.py
        - inv_helpe_functions.py
        - load_data.py
        - mesh_functions.py
        - model_selection.py
        - process_data.py
        - visualize_solution.py

### params

Contains conductivity parameters in the file **params_material.py** and specifications regarding the electrode setup in the file **params_electrode.py**.

## Important data structures

To deal with this specific problem, the following data structures were introduced:

### `params_geom`

Is a list object that contains all the geometry parameters that are necessary to create the 3D model. Its contents are:

    - params_geom
        - params_spine
        - params_vertebrae
        - params_disc
        - params_neck
        - params_trachea
        - params_oesoph

Each entry is another list that contains float values for start points, end points, radii, lenghts, etc. for the respective geometry.
The exact content is at the very end of the function `CreateParamsGeom()` inside the mesh_functions.py file.

### `data_SCP`

Is a Numpy array that contains the processed data. It is of shape (16, 1280), since 16 electrodes were used for recording and the averaged data had 1280 entries. The 1280 entries are due to the sampling rate of 512Hz and the 2.5s long signal-segment that is of interest.

### run_X.mat

Every measurement file has to contain the measurement channels as well as a channel for markers that were sent during the experiment. 
The markers used in this experiment and therefore, used in the signal analysis are:

    - 0: no marker is sent
    - 1: start run
    - 2: start trial
    - 3: start stimulation
    - 4: start post stimulation period
    - 5: end trial
    - 6: end run  

## Installations

For this software package to work, the ngsolve package has to be downloaded as well as the vtk package for visualizations.

Both packages can be downloaded here:

[Link to NGSolve](https://ngsolve.org/downloads)

[Link to VTK](https://vtk.org/download/)

A full list of all requiered packages is shown here:

| Package   |        Version       |
|-----------|:--------------------:|
|glob       |    -                 |
|math       |    -                 |
|matplotlib |    3.6.2             |
|netgen     |    6.2.2204 / 0.3.0  |
|ngsolve    |    6.2.2204          |
|numpy      |    1.23.2            |
|os         |    -                 |
|pathlib    |    1.0.1             |
|scipy      |    1.9.3             |
|sys        |    -                 |
|vtk        |    9.2.6             |


## How does it work?

### Setup

Once the data is recorded and the .mat files are as describes above, the user has to create a folder with the following name specification:

user_ + _some\_short\_ string_

The SpinalCordPotentials folder has to be located in this user_ + _some\_short\_ string_ - folder.

Additionally, the user has to create a file called measurement_geom.txt that contains the neck circumference, the Ionion-C7 distance as well as the C7-L5 distance (all in cm) with a semicolon and no whitespace in between.
For example:
- Neck circumference: 38
- Inion-C7 distance: 15
- C7-L5 distance: 50

Then the content of the measurment_geom.txt file would be:
`38;15;50`

------------

#### Working with the main.py file

Now the user can open the main.py file and set the following variables:

- `user_id`: to save the processed data, the leadfields, the mesh file, and the results in a seperate folder [string]
- `num_dipoles_per_section`: how many dipoles per interveretbral section should be assumed [int]
- `fs`: Sampling rate during measurments [int]
- `alpha`: Regularization parameter for Norm-Calculation [float]

If the user wants to perform a hyperparameter search, he can change the following boolean variables:

- `HYPERPARAMETER_SEARCH`
- `PLOT_HP_DATA`

The setup is now finished and the code can be executed.

### What happens now?

After the user presses **Run**, a 3D model is either created or loaded, depending on the geometry parameters that are given in the measurement_geom.txt file.
If no measurement_geom.txt file is provided, then the alternative geometry settings are used.

Next, the data is loaded and processed, assuming it is in the right folder (see above). 
Afterwards, the processed data is stored in the folder **processed_data**.
The user has now the possibility to plot the averaged signal including the waveform points by uncommenting the functions `ChannelPlotWP_1()` and `ChannelPlotWP_2()`.

The next step in the program is to calculate the electrode positions for the given geometry.
In this version, the electrode pattern has 16 electrodes on the right hand side of the neck.

If the parameter `HYPERPARAMETER_SEARCH` is set to **True**, the initially set variables `num_dipoles_per_section` and `alpha` are ignored and calculated according to the `HyperParameterSearch()` function.\
Furhter, the user can set the variable `PLOT_HP_DATA` to **True**, so that the Error-Functional (which is depended on the hyperparamters) is plotted.

The last step is to calcualte the inverse problem and visualize it at specific time points (P1, N1, P2).\
The intial view of the solution can be set to:

    - frontal
    - lower_frontal
    - upper_frontal
    - upper_right
    - upper_left

### Images

All figures that were created can be saved to the folder **images**, which is created automatically.

## Credits
Author: Markus E. Oberndorfer


