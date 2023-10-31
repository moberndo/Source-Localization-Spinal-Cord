""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate
    the voltage distribution caused by a dipole
    Subtitle: Visualization of the Solution
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
import vtk
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkRenderingCore import vtkWindowToImageFilter
from  math import sqrt, atan, pi
from numpy import array, zeros, mean, argmax, std, argmin
from numpy import sqrt as np_sqrt

# ----- FUNCTIONS -----
def SaveImage(filename, renWin):
    """
    This function creates the setup that is needed to save
    the vtk image.

    Input: 
        - filename: Filename under which the image should be saved. [str]
        - renWin: The renderer that holds the window which should be saved.
        [vtk render-object]
    
    Output:
        - No output.
    """
    writer = vtkPNGWriter()
    window_to_image_filter = vtkWindowToImageFilter()
    window_to_image_filter.SetInput(renWin)
    window_to_image_filter.SetScale(1)  # image quality
    window_to_image_filter.SetInputBufferTypeToRGBA()

    writer.SetFileName(filename)
    writer.SetInputConnection(window_to_image_filter.GetOutputPort())
    writer.Write()

def CreateJMagnitude(J):
    """
    From the given solution vector, every three entries represent one dipole,
    for which the total magnitude should be calcualted.

    Input:
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]

    Output:
        J_mag: Contains the magnitudes of all dipoles in J. [Numpy array]
    """
    rows, cols = J.shape[0], int(J.shape[1] / 3)
    J_mag = zeros(shape=(rows, cols))
    for i in range(cols):
        idx = i * 3
        J_mag[:,i] = np_sqrt(J[:,idx]**2 + J[:,idx+1]**2 + J[:,idx+2]**2)
    return J_mag

def GetWaveformPoints(J, data_SCP, fs):
    """
    Is a search function, that tries to find the waveform points P1, N1 and
    P2 for a given solution vector and SCP data.

    Input:
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]
        - data_SCP: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data. [Numpy array]
        - fs: Sampling frequency during experiment. [int]

    Output:
        - P1: Index of the first waveform point. [int]
        - N1: Index of the second waveform point. [int]
        - P2: Index of the third waveform point. [int]
    """

    # always plus 0.5 since the signal starts 0.5s before the stim onset
    t_idx_0 = int(fs*(0 + 0.5))
    t_idx_500 = int(fs*(0.5 + 0.5))
    t_idx_1000 = int(fs*(1 + 0.5))

    # P1
    P1_idx = [argmax(data_SCP[i, t_idx_0:t_idx_500]) + t_idx_0 \
              for i in range(16)]
    #print('P1: ',(mean(P1_idx)-t_idx_0)/fs, std(P1_idx)/fs)
    P1_ = (int(mean(P1_idx)), int(std(P1_idx)))
    
    # N1
    N1_idx = [argmin(data_SCP[i, int(mean(P1_idx)):t_idx_1000]) \
              + int(mean(P1_idx)) for i in range(16)]
    #print('N1: ',(mean(N1_idx)-t_idx_0)/fs, std(N1_idx)/fs)
    N1_ = (int(mean(N1_idx)), int(std(N1_idx)))
    
    # P2
    P2_idx = [argmax(data_SCP[i, int(mean(N1_idx)):t_idx_1000]) \
              + int(mean(N1_idx)) for i in range(16)]
    #print('P2: ',(mean(P2_idx)-t_idx_0)/fs, std(P2_idx)/fs)
    P2_ = (int(mean(P2_idx)), int(std(P2_idx)))

    J_mag = CreateJMagnitude(J)

    P1 = [argmax(J_mag[P1_[0] - P1_[1]:P1_[0] + P1_[1],i]) \
          + P1_[0] - P1_[1] for i in range(J_mag.shape[1])]
    P1 = int(mean(P1))
    N1 = [argmax(J_mag[N1_[0] - N1_[1]:N1_[0] + N1_[1],i]) \
          + N1_[0] - N1_[1] for i in range(J_mag.shape[1])]
    N1 = int(mean(N1))
    P2 = [argmax(J_mag[P2_[0] - P2_[1]:P2_[0] + P2_[1],i])  \
          + P2_[0] - P2_[1] for i in range(J_mag.shape[1])]
    P2 = int(mean(P2))
    
    return P1, N1, P2

def ShowSolution(J, params_geom, z_vals, num_dipoles_per_section,
                 view = 'frontal', save = False, path = None):
    ''' 
     Visualization of the solution vector J using Visualization
     Toolkit VTK.

     Input:
        - J: Solution vector that contains the magnitude for every dipole. 
        [list, float]
        - params_geom: software intern structure. See ReadMe.txt file
          for more information. [-]
        - z_vals: z values for all dipoles [list; int]
        - num_dipoles_per_section: number of dipoles per intervertebral
          section that are assumed a priori. [int] 
        - view: Settings that defines the intial view onto the solution. 
        Possible values: 'frontal', 'lower_frontal', 'upper_frontal',
        'upper_right', 'upper_left'. [str]
        - save: Boolean variable that defines wether the solution should
        be saved. [bool]
        - path: If 'save' is set to True, the path has to be set. [str]

    Output:
        - No output.
    '''
    def CreateRenderObjects():
        """
        Set up the intial empty render object.

        Input:
            - No input.
        Output:
            - ren: Renderer object [VTK renderer object]
            - renWin: RenderWindow in which the solution is shown.
            [VTK renderwindow object]
            - iren: Interactor that is needed to manipulate the objects
             inside the RenderWindow. [VTK Interactor object]
        """

        # Create a rendering window and renderer
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)

        # Set the desired window size (width, height)
        window_width = 1800
        window_height = 800
        renWin.SetSize(window_width, window_height)

        # Create an interactor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        style = vtk.vtkInteractorStyleTrackballCamera()
        iren.SetInteractorStyle(style)

        return ren, renWin, iren
    
    def CreateArrowData(J, z_vals, neck_rad, params_geom, y_shift):
        """
        For the given solution vector and geometry, this function
        calculates the dipoles vectors.

        Input:

            - J: Solution vector that contains the magnitude for every dipole.
            [list, float]
            - z_vals: z values for all dipoles [list; int]
            - neck_rad: Radius of the neck used for this geoemtry. [float]
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - arrow_start_points: Array that contain all the start points of
            the dipole arrows. [Numpy array]
            arrow_end_points: Array that contain all the end points of
            the dipole arrows. [Numpy array]
        """
        #rad_dipole = params_geom[0][1] * 0.25
        arrow_start_points = array([[0, y_shift, z_vals[i]] \
                                    for i in range(z_vals.shape[0])])
        arrow_direction = [[J[int(i*3)], J[int(i*3+1)], J[int(i*3+2)]] \
                           for i in range(int(J.shape[0]/3))]
        
        # Scale arrows, such that the maximal length is not
        # larger then neck_rad
        arrow_lengths = []
        max_arrow_length = 0
        for _, arrow in enumerate(arrow_direction):
            arrow_len = vec_len(arrow)
            if arrow_len > max_arrow_length:
                max_arrow_length = arrow_len
                scaling_factor = neck_rad / arrow_len
            arrow_lengths.append(arrow_len)
        arrow_direction = array(arrow_direction) * scaling_factor
        
        arrow_end_points = arrow_start_points + arrow_direction
        return arrow_start_points, arrow_end_points

    def CreateSpinalCordActor(params_geom, y_shift):
        """
        Creates an vtk-actor for the spinal cord.

        Input:
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - spinal_cord_actor: Cylinder object that represents
            the spinal cord in the solution visualization. 
            [VTK actor object]
        """
        # Spinal Cord Params
        spinal_cord_height = params_geom[3][2]
        spinal_cord_rad = params_geom[0][1]
        spinal_cord_start = params_geom[0][0]
        spinal_cord_end = spinal_cord_start - spinal_cord_height

        # Create Neck
        spinal_cord = vtk.vtkCylinderSource()
        spinal_cord.SetResolution(100)
        spinal_cord.SetHeight(spinal_cord_height)
        spinal_cord.SetRadius(spinal_cord_rad)
        spinal_cord.Update()

        # Create a mapper and actor for the spinal_cord
        spinal_cord_mapper = vtk.vtkPolyDataMapper()
        spinal_cord_mapper.SetInputData(spinal_cord.GetOutput())
        spinal_cord_actor = vtk.vtkActor()
        spinal_cord_actor.GetProperty().SetOpacity(0.4)
        spinal_cord_actor.GetProperty().SetColor(0.55, 0.45, 0.45)
        spinal_cord_actor.SetMapper(spinal_cord_mapper)

        # Transform the cylinder
        spinal_cord_transform = vtk.vtkTransform()
        spinal_cord_transform.Translate([0,
                                         y_shift,
                                         (spinal_cord_end \
                                          + spinal_cord_height / 2)])
        spinal_cord_transform.RotateX(-90.)
        spinal_cord_actor.SetUserTransform(spinal_cord_transform)

        return spinal_cord_actor
    
    def CreateNeckActor(params_geom, y_shift):
        """
        Creates an vtk-actor for the neck.

        Input:
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - neck_actor: Cylinder object that represents the neck in the
            solution visualization. [VTK actor object]
        """
        # Neck Params
        neck_dist = params_geom[3][1]
        neck_height = params_geom[3][2]
        neck_rad = params_geom[3][0]
        neck_start = params_geom[0][0]
        neck_end = neck_start - neck_height

        # Create Neck
        neck = vtk.vtkCylinderSource()
        neck.SetResolution(100)
        neck.SetHeight(neck_height)
        neck.SetRadius(neck_rad)
        neck.Update()

        # Create a mapper and actor for the neck
        neck_mapper = vtk.vtkPolyDataMapper()
        neck_mapper.SetInputData(neck.GetOutput())
        neck_actor = vtk.vtkActor()
        neck_actor.GetProperty().SetOpacity(0.3)
        neck_actor.SetMapper(neck_mapper)

        # Transform the cylinder
        neck_transform = vtk.vtkTransform()
        neck_transform.Translate([neck_dist,
                                  y_shift,
                                  (neck_end + neck_height / 2)])
        neck_transform.RotateX(-90.)
        neck_actor.SetUserTransform(neck_transform)

        return neck_actor

    def CreateArrowActors(J, z_vals, params_geom, y_shift):
        """
        Creates an vtk-actor for the solution arrows.

        Input:
            - J: Solution vector that contains the magnitude for every dipole.
            [list, float]
            - z_vals: z values for all dipoles [list; int]
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - arrow_actors: List of arrow objects that represents the arrows
            in the solution visualization. [VTK actor object]
        """
        # Create arrows
        arrow_source = vtk.vtkArrowSource()
        arrow_source_dir = [1,0,0]

        # Define arrow start and end points
        arrow_out_ = CreateArrowData(J, z_vals, params_geom[3][0], params_geom, y_shift)
        arrow_start_points = arrow_out_[0]
        arrow_end_points = arrow_out_[1]
        # Create a list to hold the arrow actors
        arrow_actors = []
        
        # Create actors for arrows
        for start, end in zip(arrow_start_points, arrow_end_points):
            arrow = vtk.vtkActor()
            arrow_mapper = vtk.vtkPolyDataMapper()
            arrow_mapper.SetInputConnection(arrow_source.GetOutputPort())
            arrow.SetMapper(arrow_mapper)
            arrow.GetProperty().SetColor(1.0, 0.0, 0.0)
    
            # Calculate arrow transformation
            arrow_transform = vtk.vtkTransform()
            arrow_transform.Translate(start)

            arrow_direction = [end[i] - start[i] for i in range(3)]
            arrow_len = vec_len(arrow_direction)

            arrow_z_len = arrow_direction[2]

            arrow_xy = [arrow_direction[0],arrow_direction[1], 0]
            arrow_xy_len = vec_len(arrow_xy)
            
            alpha = atan(arrow_z_len / arrow_xy_len) * 180 / pi

            if arrow_xy_len == 0.:
                alpha = 90
                beta = 0
            # check for first quadrant
            elif arrow_direction[0] > 0 and arrow_direction[1] > 0: 
                x_comp = arrow_direction[0]
                y_comp = arrow_direction[1]
                beta = atan(y_comp / x_comp) * 180 / pi
            # check for second quadrant
            elif arrow_direction[0] < 0 and arrow_direction[1] > 0:      
                x_comp = -arrow_direction[0]
                y_comp = arrow_direction[1]
                beta = 180 - atan(y_comp / x_comp) * 180 / pi
            # check for third quadrant
            elif arrow_direction[0] < 0 and arrow_direction[1] < 0:      
                x_comp = -arrow_direction[0]
                y_comp = -arrow_direction[1]
                beta = 180 + atan(y_comp / x_comp) * 180 / pi
            # check for fourth quadrant
            else:                                                        
                x_comp = arrow_direction[0]
                y_comp = -arrow_direction[1]
                beta = 360 - atan(y_comp / x_comp) * 180 / pi

            arrow_transform.RotateY(alpha)
            arrow_transform.RotateZ(beta)
            arrow_transform.Scale(arrow_len,arrow_len,arrow_len)

            arrow.SetUserTransform(arrow_transform)
            arrow_actors.append(arrow)
        return arrow_actors
    
    def AddTitle(ren, y_shift, title, params_geom):
        """
        Adds a titel to every cylinder that represents the neck geometry
        at a differnt time point (P1, N1, P2).

        Input:
            - ren: Renderer object [VTK renderer object]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]
            - title: Title that should be applied. [str]
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]

        Output:
            - No output.
        """
        # define params
        neck_height = params_geom[3][2]
        neck_start = params_geom[0][0]
        spinal_cord_rad = params_geom[0][1]
        loc = [0, y_shift - spinal_cord_rad*1.25,
               neck_start + neck_height*0.2]

        #create actor
        title_text = vtk.vtkVectorText()
        title_text.SetText(title)

        title_text_mapper = vtk.vtkPolyDataMapper()
        title_text_mapper.SetInputConnection(title_text.GetOutputPort())

        title_text_actor = vtk.vtkActor()
        title_text_actor.SetMapper(title_text_mapper)
        title_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)

        title_text_transform = vtk.vtkTransform()
        title_text_transform.Translate(loc)
        title_text_transform.Scale(0.0175, 0.0175, 0.0175)
        title_text_transform.RotateZ(90)
        title_text_transform.RotateX(90)
        
        title_text_actor.SetUserTransform(title_text_transform)

        # add actor to render window
        ren.AddActor(title_text_actor)

    def CreateCoordinateArrows(params_geom, y_shift):
        """
        Create a coordinate system that is placed beside the three 
        visualizations of the solution.

        Input: 
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - arrow_actors: a list of all actor objects that represent
            the arrows. [list, VTK actor object]
            - arrow_text_actors: A list of all texts that are shown beside
            the coordinate arrows. [list, VTK actor object]
        """
        # Params
        neck_height = params_geom[3][2]
        neck_rad = params_geom[3][0]
        neck_start = params_geom[0][0]
        neck_end = neck_start - neck_height
        arrows_len = neck_rad

        start = [0, y_shift, neck_end] # (neck_end + neck_height / 2)
        end_right = [0, y_shift + arrows_len*0.9, neck_end*1.025]
        #end_left = [0, y_shift-arrows_len*1.65, neck_end*1.025]
        end_forward = [arrows_len, y_shift*1.02, neck_end]
        end_upper = [0, y_shift-arrows_len*0.25, arrows_len*1.2 + neck_end]

        # Create arrows
        arrow_source = vtk.vtkArrowSource()
        arrow_actors = []
        arrow_text_actors = []
        
        # Create actors for arrows
        for _ in range(3):
            arrow = vtk.vtkActor()
            arrow_mapper = vtk.vtkPolyDataMapper()
            arrow_mapper.SetInputConnection(arrow_source.GetOutputPort())
            arrow.SetMapper(arrow_mapper)
            arrow.GetProperty().SetColor(0.15, 0.15, 0.15)
            arrow_actors.append(arrow)

        '''
        # right - arrow transformation
        right_arrow_transform = vtk.vtkTransform()
        right_arrow_transform.Translate(start)
        #right_arrow_transform.Scale(arrows_len, arrows_len, arrows_len)
        right_arrow_transform.RotateZ(90)
        arrow_actors[0].SetUserTransform(right_arrow_transform)

        # right - arrow text
        right_arrow_text = vtk.vtkTextActor()
        right_arrow_text.SetInput("Dexter")

        right_arrow_text_prop = vtk.vtkTextProperty()
        right_arrow_text_prop.SetFontFamilyToTimes()
        right_arrow_text_prop.SetFontSize(10)
        #right_arrow_text.SetItalic(True)

        right_arrow_text.SetTextProperty(right_arrow_text_prop)

        right_arrow_text_mapper = vtk.vtkPolyDataMapper()
        right_arrow_text_mapper.SetInputConnection(right_arrow_text.GetOutputPort())

        #right_arrow_text_actor = vtk.vtkActor()
        #right_arrow_text_actor.GetProperty().SetTextProperty()
        #right_arrow_text_actor.SetTextProperty(right_arrow_text_prop)
        #right_arrow_text_actor.SetMapper(right_arrow_text_mapper)
        #right_arrow_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)
        #right_arrow_text_actor.GetProperty().SetItalic(True)

        right_arrow_text_transform = vtk.vtkTransform()
        right_arrow_text_transform.Translate(end_right)
        right_arrow_text_transform.Scale(0.0175, 0.0175, 0.0175)
        right_arrow_text_transform.RotateZ(90)
        right_arrow_text_transform.RotateX(90)
        
        right_arrow_text.SetUserTransform(right_arrow_text_transform)
        

        arrow_text_actors.append(right_arrow_text)
        '''
        # left - arrow transformation
        left_arrow_transform = vtk.vtkTransform()
        left_arrow_transform.Translate(start)
        left_arrow_transform.Scale(arrows_len, arrows_len, arrows_len)
        left_arrow_transform.RotateZ(90)
        arrow_actors[0].SetUserTransform(left_arrow_transform)

        # left - arrow text
        left_arrow_text = vtk.vtkVectorText()
        left_arrow_text.SetText("Dexter")

        left_arrow_text_mapper = vtk.vtkPolyDataMapper()
        left_arrow_text_mapper.SetInputConnection(left_arrow_text.GetOutputPort())
        
        left_arrow_text_actor = vtk.vtkActor()
        left_arrow_text_actor.SetMapper(left_arrow_text_mapper)
        left_arrow_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)

        left_arrow_text_transform = vtk.vtkTransform()
        left_arrow_text_transform.Translate(end_right)
        left_arrow_text_transform.Scale(0.0175, 0.0175, 0.0175)
        left_arrow_text_transform.RotateZ(90)
        left_arrow_text_transform.RotateX(90)
        
        left_arrow_text_actor.SetUserTransform(left_arrow_text_transform)

        arrow_text_actors.append(left_arrow_text_actor)
        

        '''
        # left - arrow transformation
        left_arrow_transform = vtk.vtkTransform()
        left_arrow_transform.Translate(start)
        left_arrow_transform.Scale(arrows_len, arrows_len, arrows_len)
        left_arrow_transform.RotateZ(-90)
        arrow_actors[1].SetUserTransform(left_arrow_transform)

        # left - arrow text
        left_arrow_text = vtk.vtkVectorText()
        left_arrow_text.SetText("Left")

        left_arrow_text_mapper = vtk.vtkPolyDataMapper()
        left_arrow_text_mapper.SetInputConnection(left_arrow_text.GetOutputPort())
        
        left_arrow_text_actor = vtk.vtkActor()
        left_arrow_text_actor.SetMapper(left_arrow_text_mapper)
        left_arrow_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)

        left_arrow_text_transform = vtk.vtkTransform()
        left_arrow_text_transform.Translate(end_left)
        left_arrow_text_transform.Scale(0.0175, 0.0175, 0.0175)
        left_arrow_text_transform.RotateZ(90)
        left_arrow_text_transform.RotateX(90)
        
        left_arrow_text_actor.SetUserTransform(left_arrow_text_transform)

        arrow_text_actors.append(left_arrow_text_actor)
        '''

        # forward - arrow transformation
        forward_arrow_transform = vtk.vtkTransform()
        forward_arrow_transform.Translate(start)
        forward_arrow_transform.Scale(arrows_len, arrows_len, arrows_len)
        #forward_arrow_transform.RotateZ(-90)
        arrow_actors[1].SetUserTransform(forward_arrow_transform)

        # forward - arrow text
        forward_arrow_text = vtk.vtkVectorText()
        forward_arrow_text.SetText("Anterior")

        forward_arrow_text_mapper = vtk.vtkPolyDataMapper()
        forward_arrow_text_mapper.SetInputConnection(forward_arrow_text.GetOutputPort())
        
        forward_arrow_text_actor = vtk.vtkActor()
        forward_arrow_text_actor.SetMapper(forward_arrow_text_mapper)
        forward_arrow_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)

        forward_arrow_text_transform = vtk.vtkTransform()
        forward_arrow_text_transform.Translate(end_forward)
        forward_arrow_text_transform.Scale(0.0175, 0.0175, 0.0175)
        forward_arrow_text_transform.RotateZ(90)
        forward_arrow_text_transform.RotateX(90)
        
        forward_arrow_text_actor.SetUserTransform(forward_arrow_text_transform)

        arrow_text_actors.append(forward_arrow_text_actor)

        # upper - arrow transformation
        upper_arrow_transform = vtk.vtkTransform()
        upper_arrow_transform.Translate(start)
        upper_arrow_transform.Scale(arrows_len, arrows_len, arrows_len)
        upper_arrow_transform.RotateZ(90)
        upper_arrow_transform.RotateY(-90)
        arrow_actors[2].SetUserTransform(upper_arrow_transform)

        # upper - arrow text
        upper_arrow_text = vtk.vtkVectorText()
        upper_arrow_text.SetText("Cranial")

        upper_arrow_text_mapper = vtk.vtkPolyDataMapper()
        upper_arrow_text_mapper.SetInputConnection(upper_arrow_text.GetOutputPort())
        
        upper_arrow_text_actor = vtk.vtkActor()
        upper_arrow_text_actor.SetMapper(upper_arrow_text_mapper)
        upper_arrow_text_actor.GetProperty().SetColor(0.15, 0.15, 0.15)

        upper_arrow_text_transform = vtk.vtkTransform()
        upper_arrow_text_transform.Translate(end_upper)
        upper_arrow_text_transform.Scale(0.0175, 0.0175, 0.0175)
        upper_arrow_text_transform.RotateZ(90)
        upper_arrow_text_transform.RotateX(90)
        
        upper_arrow_text_actor.SetUserTransform(upper_arrow_text_transform)

        arrow_text_actors.append(upper_arrow_text_actor)

        return arrow_actors, arrow_text_actors

    def CreateSolution(params_geom, J, z_vals, y_shift):
        """
        Creates all actors that are necessary for the visualization.

        Input:
            - params_geom: software intern structure. See ReadMe.txt file
            for more information. [-]
            - J: Solution vector that contains the magnitude for every dipole.
            [list, float]
            - z_vals: z values for all dipoles [list; int]
            - y_shift: Shift between the Neck-Cylinder and 
            the SpinalCord-Cylinder. [float]

        Output:
            - neck_actor: Actor for the neck visualiuation. [VTK actor object]
            - spinal_cord_actor: Actor for the spinal cord visualiuation.
            [VTK actor object]
            - arrow_actors: Actor for the solution vector visualiuation.
            [VTK actor object]
        """
        neck_actor = CreateNeckActor(params_geom, y_shift)
        spinal_cord_actor = CreateSpinalCordActor(params_geom, y_shift)
        arrow_actors = CreateArrowActors(J, z_vals, params_geom, y_shift)

        return neck_actor, spinal_cord_actor, arrow_actors

    def AddSolutionToRen(ren, actors):
        """
        Takes the actors and adds them to a given VTK renderer.

        Input:
            - ren: Renderer object [VTK renderer object]
            - actors: List of all actors that should be added.
            [list, VTK actor object]

        Output:
            - No output.
        """
        
        neck_actor = actors[0]
        spinal_cord_actor = actors[1]
        arrow_actors = actors[2]

        ren.AddActor(neck_actor)
        ren.AddActor(spinal_cord_actor)
        for arrow_actor in arrow_actors:
            ren.AddActor(arrow_actor)

    # define parameter
    neck_rad = params_geom[3][0]
    y_shift = neck_rad * 3
    # Calculate vector length
    vec_len = lambda vec: sqrt( sum( vec[i]**2 for i in range(3) ) )

    ''' MAIN '''
    ren, renWin, iren = CreateRenderObjects()
    # solution for P1
    J_P1 = J[0]
    y_shift_P1 = y_shift * 0
    actors_P1 = CreateSolution(params_geom, J_P1, z_vals, y_shift_P1)
    AddSolutionToRen(ren, actors_P1)
    AddTitle(ren, y_shift_P1, 'P1', params_geom)

    # solution for N1
    J_N1 = J[1]
    y_shift_N1 = y_shift * 1
    actors_N1 = CreateSolution(params_geom, J_N1, z_vals, y_shift_N1)
    AddSolutionToRen(ren, actors_N1)
    AddTitle(ren, y_shift_N1, 'N1', params_geom)

    # solution for P2
    J_P2 = J[2]
    y_shift_P2 = y_shift * 2
    actors_P2 = CreateSolution(params_geom, J_P2, z_vals, y_shift_P2)
    AddSolutionToRen(ren, actors_P2)
    AddTitle(ren, y_shift_P2, 'P2', params_geom)

    # Add Coordinate System
    y_shift_CS = y_shift * 2.5
    out_ = CreateCoordinateArrows(params_geom, y_shift_CS)
    coordinate_arrow_actors = out_[0]
    coordinate_arrow_texts = out_[1]
    for coordinate_arrow_actor in coordinate_arrow_actors:
            ren.AddActor(coordinate_arrow_actor)
    for coordinate_arrow_text in coordinate_arrow_texts:
        ren.AddActor(coordinate_arrow_text)

    # Set the background color
    ren.SetBackground(0.1, 0.1, 0.1)

    # Set up the camera and renderer
    camera_transform = vtk.vtkTransform()

    # change view
    camera_transform.RotateWXYZ(90, 0, 1, 0)
    camera_transform.RotateWXYZ(90, 0, 0, 1)
    camera_transform.RotateWXYZ(-20, 1, 0, 0)

    if view == 'frontal':
        filename = 'NV_' + str(num_dipoles_per_section)
        filename = filename + '_frontal.png'
    elif view == 'upper_left':
        filename = 'NV_' + str(num_dipoles_per_section)
        filename = filename + '_upper_left.png'
        camera_transform.RotateWXYZ(-25, 0, 1, 0)
        camera_transform.RotateWXYZ(-5, 1, 0, 0)
    elif view == 'upper_right':
        filename = 'NV_' + str(num_dipoles_per_section)
        filename = filename + '_upper_right.png'
        camera_transform.RotateWXYZ(25, 0, 1, 0)
        camera_transform.RotateWXYZ(-5, 1, 0, 0)
    elif view == 'upper_frontal':
        filename = 'NV_' + str(num_dipoles_per_section)
        filename = filename + '_upper_frontal.png'
        camera_transform.RotateWXYZ(-55, 1, 0, 0)
    elif view == 'lower_frontal':
        filename = 'NV_' + str(num_dipoles_per_section)
        filename = filename + '_lower_frontal.png'
        camera_transform.RotateWXYZ(55, 1, 0, 0)
        

    ren.GetActiveCamera().ApplyTransform(camera_transform)
    ren.GetActiveCamera().Azimuth(0)
    ren.GetActiveCamera().Elevation(0)
    
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(2.2)
    # Start the interaction
    iren.Initialize()
    renWin.Render()
    iren.Start()

    #return ren, iren, renWin
