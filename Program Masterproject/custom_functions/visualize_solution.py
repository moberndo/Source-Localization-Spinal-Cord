import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.animation import FuncAnimation
from numpy import array, linspace, amax, mean
from functools import partial
#from params import CreateParamsGeom

def GetPntsPoly(y0, height_vert, rad_spine):

    pnts_poly = array([[-rad_spine*1.1, -rad_spine*1.1, -rad_spine*1.1 - 0.0125, -rad_spine*1.1 - 0.015],
                       [y0+height_vert*0.4, y0+height_vert*0.9, y0+height_vert*0.6, y0+height_vert*0.2]])

    return pnts_poly

def DrawBackground(params_geom, z_vals):

    loc_neck = params_geom[3][1]
    rad_neck = params_geom[3][0]
    rad_spine = params_geom[0][1]
    len_spine = params_geom[0][0]
    rad_vertebrae = params_geom[1][0]
    height_vertebrae = params_geom[1][3]
    height_vert_thoracic = params_geom[1][5]
    height_disc =params_geom[2][0]
    
    
    fig, ax = plt.subplots()

    # add neck outlines
    ax.plot([loc_neck-rad_neck, loc_neck-rad_neck],
            [len_spine-7*(height_disc+height_vertebrae)-height_vert_thoracic, len_spine], color = 'black', linewidth = 3.5)
    ax.plot([loc_neck+rad_neck, loc_neck+rad_neck],
            [len_spine-7*(height_disc+height_vertebrae)-height_vert_thoracic, len_spine], color = 'black', linewidth = 3.5)

    # add spine outlines
    ax.plot([-rad_spine, -rad_spine],
            [len_spine-7*(height_disc+height_vertebrae)-height_vert_thoracic, len_spine], color = 'black', linewidth = 2)
    ax.plot([rad_spine, rad_spine],
            [len_spine-7*(height_disc+height_vertebrae)-height_vert_thoracic, len_spine], color = 'black', linewidth = 2)

    # add vertebraes
    y0_OG = len_spine - (height_disc + height_vertebrae)
    for vert_idx in range(7):   # change to 6 to have no disc between atlas and axis
        y0 = y0_OG - vert_idx * (height_disc + height_vertebrae)
        ax.add_patch(Rectangle( (rad_spine*1.1, y0), width = rad_vertebrae*1.85, height = height_disc,
                     edgecolor = 'gray', facecolor = 'gray', fill = True))
        
    for vert_idx in range(7):
        '''
        ##########    remove comments to have no disc between atlas and axis      ############
        if vert_idx == 6:
            y0 = y0_OG - vert_idx * (height_disc + height_vertebrae)
        else:
            y0 = y0_OG - vert_idx * (height_disc + height_vertebrae) + height_disc'''
        
        y0 = y0_OG - vert_idx * (height_disc + height_vertebrae) + height_disc

        ax.add_patch(Rectangle( (rad_spine*1.1, y0), width = rad_vertebrae*1.85, height = height_vertebrae,
                     edgecolor = 'black', fill = False, lw = 1.5))
        C_idx = 'C'+str(vert_idx+1)
        ax.text(rad_spine*1.1+rad_vertebrae, y0+height_vertebrae*0.5, C_idx, style = 'oblique', 
                verticalalignment='center', horizontalalignment='right')
        
        pnts_poly = GetPntsPoly(y0, height_vertebrae, rad_spine)
        ax.add_patch(Polygon(xy = pnts_poly.T, closed=True, edgecolor = 'black', fill = False, lw = 1.5))

    # add thoracic vertebrae
    y0 = y0_OG - 6 * (height_disc + height_vertebrae) - height_vert_thoracic
    ax.add_patch(Rectangle( (rad_spine*1.1, y0), width = rad_vertebrae*1.85, height = height_vert_thoracic,
                     edgecolor = 'black', fill = False, lw = 1.5 ))
    ax.text(rad_spine*1.1+rad_vertebrae, y0+height_vert_thoracic*0.5, 'Th1', style = 'oblique', 
            verticalalignment='center', horizontalalignment='right')
    
    pnts_poly = GetPntsPoly(y0, height_vert_thoracic, rad_spine)
    ax.add_patch(Polygon(xy = pnts_poly.T, closed=True, edgecolor = 'black', fill = False, lw = 1.5))

    # add solution elements
    if len(z_vals) >= 14:
        height_uni = abs(z_vals[1]- z_vals[0]) / 2
    else:
        height_uni = height_vertebrae * 0.25

    for _, val in enumerate(z_vals):
        ax.add_patch(Rectangle( (-rad_spine, val - height_uni), width = rad_spine*2, height = height_uni*2,
                            edgecolor = 'white', facecolor = 'white', fill = True, alpha = 0, lw = 0,
                            label = str(val)))

    plt.title('Dipoles Right')

    patches_ = [patch for patch in ax.patches if patch.get_label() != '']
    #plt.axis('off')

    return ax, fig, patches_

def DrawSolution(i, J, z_vals):

    J_scale = J[i,:,0] / max_val
    '''
    vals_per_section = int(len(z_vals) / 7)

    for j in range(7):
        mean_val = mean(J_scale[vals_per_section * j : vals_per_section * j + vals_per_section])
        J_scale[vals_per_section * j : vals_per_section * j + vals_per_section] = mean_val
        
    '''

    for idx, val in enumerate(z_vals):
        for patch in patches:
            if patch.get_label() == str(val):
                if J_scale[idx] >= 0:
                    patch.set_facecolor('red')
                    patch.set_edgecolor('red')
                    patch.set_alpha(J_scale[idx])
                elif J_scale[idx] < 0:
                    patch.set_facecolor('blue')
                    patch.set_edgecolor('blue')
                    patch.set_alpha(-J_scale[idx])

    return patches

def init():
    for patch in patches:
        patch.set_facecolor('white')
        patch.set_edgecolor('white')
    return patches

def CreateAnimation(J, params_geom, z_vals):
    _, fig_back, patches_ = DrawBackground(params_geom, z_vals)

    global max_val
    max_val = amax(abs(J))
    global patches 
    patches = patches_
    n_frames = J.shape[0]
    #print(n_frames)

    animator = FuncAnimation(fig_back, partial(DrawSolution, J = J, z_vals = z_vals),
                             init_func= init,
                             frames = range(n_frames),
                             interval = 1 / 0.512)

    return animator


if __name__ == '__main__':
    rad_neck = 0.06
    loc_neck = 0.03
    rad_spine = 0.01
    loc_C7 = 0.54
    len_spine = 0.675
    rad_vertebrae = 0.0175
    height_vertebrae = 0.015
    height_disc = 0.003
    height_vert_thoracic = 0.022

    #params_geom = CreateParamsGeom('39;15;50')

    J = array([-1,0,0,0,0,0,0.1,0.2,0.4,0.8,0,0,0,1,1.4,1.7,4,3.6,2.1,0.6,0.1,0,0,0,0,0.1,0.1,0.1,0.2,0.4,0,0,0.1,-1])

    #ax, fig = DrawBackground(params_geom)
    #DrawSolution(ax, J, 1)