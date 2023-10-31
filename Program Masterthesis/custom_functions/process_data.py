""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate
    the voltage distribution caused by a dipole
    Subtitle: Process data Function
    Author: Markus E. Oberndorfer
####################################################################### """

# ----- IMPORT PACKAGES -----
from scipy.fft import fft
from scipy.signal import filtfilt, butter, iirnotch
from numpy import ones_like, array, linspace, where, zeros, mean, std
from numpy import max, min, save, argmax, argmin
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec

# ----- FUNCTIONS -----
def ChannelPlotWP_1(data, fs, save = False, main_path = None, user_ID = None):
    """
    Creates the Channel plot that visualizes the averaged data for one
    electrode positions and marks the waveform points.

    Input:
        - data: Averaged SCP data. [Numpy array]
        - fs: Sampling frequency. [int]
        - save: If set to True, the image is saved to the given user_ID path.
        [bool]
        - main_path: Contains the path of the main file [pathlib object]
        - user_ID: Is a short ID that points to the folder structure of
        interest. [str]
    
    Output:
        - No output.
    """

    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams['lines.linewidth'] = 2.0

    _, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (6,5))
    ax.set_yticks([-2, 0, 2], ['-2', '0', '2'])
    ax.set_ylabel('amplitude in $\mu$V', fontsize = 15)
    ax.set_xticks([0, 1, 2],['0', '1', '2'])
    ax.set_xlabel('time in s', fontsize = 15)
    ax.set_title('Sp1', fontweight = 'bold', fontsize = 15)
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis = 'both', which = 'both', length = 0)
    ax.set_ylim([-2.5, 2.5])

    # define x-axis data
    t = linspace(-0.5, 2, int(fs*2.5))
    # define y-axis data
    y = data[5,:]

    # find P1, N1 and P2
    P1 = argmax(y[0:fs])
    P1 = (P1, y[P1])
    N1 = argmin(y[P1[0]:int(fs*2)]) + P1[0]
    N1 = (N1, y[N1])
    P2 = argmax(y[N1[0]:]) + N1[0]
    P2 = (P2, y[P2])
    

    # plot data and helper lines
    ax.plot(t, y, 'k', alpha = 1, linewidth = '0.8')
    ax.vlines(0, -3, 3, colors = 'grey', linestyles = 'dotted')
    ax.vlines(1, -3, 3, colors = 'grey', linestyles = 'dotted')

    # plot waveform points
    scatter_x = [P1[0]/fs - 0.5, N1[0]/fs - 0.5, P2[0]/fs - 0.5]
    scatter_y = [P1[1], N1[1], P2[1]]
    ax.scatter(scatter_x, scatter_y, c='red', marker='x')

    ax.text(scatter_x[0] - 0.055, scatter_y[0] + 0.1, 'P1',
            style='italic', color = 'red')
    ax.text(scatter_x[1] - 0.07, scatter_y[1] - 0.25, 'N1',
            style='italic', color = 'red')
    ax.text(scatter_x[2] - 0.055, scatter_y[2] + 0.1, 'P2',
            style='italic', color = 'red')

    if save:
        user_folder = 'user_' + user_ID
        path_fig = main_path / 'data' / user_folder / 'images'
        path_fig = path_fig / 'WP_1.png'
        plt.savefig(path_fig)

def ChannelPlotWP_2(data, fs, save = False, main_path = None, user_ID = None):
    """
    Creates the Channel plot that visualizes the averaged data for one
    electrode positions and marks the waveform points.

    Input:
        - data: Averaged SCP data. [Numpy array]
        - fs: Sampling frequency. [int]
        - save: If set to True, the image is saved to the given user_ID path.
        [bool]
        - main_path: Contains the path of the main file [pathlib object]
        - user_ID: Is a short ID that points to the folder structure of
        interest. [str]
    
    Output:
        - No output.
    """

    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams['lines.linewidth'] = 2.0

    _, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (6,5))
    ax.set_yticks([-2, 0, 2], ['-2', '0', '2'])
    ax.set_ylabel('amplitude in $\mu$V', fontsize = 15)
    ax.set_xticks([0, 1, 2],['0', '1', '2'])
    ax.set_xlabel('time in s', fontsize = 15)
    ax.set_title('Sp16', fontweight = 'bold', fontsize = 15)
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(axis = 'both', which = 'both', length = 0)
    ax.set_ylim([-2.5, 2.5])

    # define x-axis data
    t = linspace(-0.5, 2, int(fs*2.5))
    # define y-axis data
    y = data[10,:]

    # find N1, P1, N2 and P2
    P1 = argmax(y[0:fs])
    P1 = (P1, y[P1])
    N1 = argmin(y[0:fs])
    N1 = (N1, y[N1])
    N2 = argmin(y[P1[0]:int(fs*2)]) + P1[0]
    N2 = (N2, y[N2])
    P2 = argmax(y[N2[0]:]) + N2[0]
    P2 = (P2, y[P2])
    

    # plot data and helper lines
    ax.plot(t, y, 'k', alpha = 1, linewidth = '0.8')
    ax.vlines(0, -3, 3, colors = 'grey', linestyles = 'dotted')
    ax.vlines(1, -3, 3, colors = 'grey', linestyles = 'dotted')

    # plot waveform points
    scatter_x = [N1[0]/fs - 0.5, P1[0]/fs - 0.5, N2[0]/fs - 0.5,
                 P2[0]/fs - 0.5]
    scatter_y = [N1[1], P1[1], N2[1], P2[1]]
    ax.scatter(scatter_x, scatter_y, c='red', marker='x')

    ax.text(scatter_x[0] - 0.07, scatter_y[0] - 0.25, 'N1$^*$',
            style='italic', color = 'red')
    ax.text(scatter_x[1] - 0.055, scatter_y[1] + 0.1, 'P1$^*$',
            style='italic', color = 'red')
    ax.text(scatter_x[2] - 0.07, scatter_y[2] - 0.25, 'N2$^*$',
            style='italic', color = 'red')
    ax.text(scatter_x[3] - 0.055, scatter_y[3] + 0.1, 'P2$^*$',
            style='italic', color = 'red')

    if save:
        user_folder = 'user_' + user_ID
        path_fig = main_path / 'data' / user_folder / 'images'
        path_fig = path_fig / 'WP_2.png'
        plt.savefig(path_fig)

def Create16ChannelPlot(data, data_runs, fs):
    """
    Creates a plot that shows the averaged data for all electrodes.
    
    Input: 
        - data: Array of shape (NumberElectrodes x SignalLength) that
        contains the processed data, that is averages over all
        participants. [Numpy array]
        - data_runs: Array that contains the data for every particiapant.
        [Numpy array]
        - fs: Sampling rate during the measurment. [int]

    Output:
        - No Output.
    """
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12

    # general figure settings
    fig, ax_ = plt.subplots(ncols=1, nrows=1, figsize = (10,12))
    ax_.set_yticks([], [])
    ax_.set_ylabel('amplitude in $\mu$V', fontsize = 15)
    ax_.set_xticks([],[])
    ax_.set_xlabel('time in s', fontsize = 15)
    ax_.spines[['right', 'top']].set_visible(False)
    ax_.spines['left'].set_position(('axes', -0.05))
    ax_.spines['left'].set_bounds(-0.05, 1)
    ax_.spines['bottom'].set_position(('axes', -0.05))
    ax_.spines['bottom'].set_bounds(-0.05, 1)
    ax_.text(-0.0925, 0.165, 'C7', fontweight = 'bold', fontsize = 15)
    ax_.text(-0.0875, 0.145, r'$\mathbf{\times}$',
             fontweight = 'bold', fontsize = 18)
    
    gs = GridSpec(12, 3, figure = fig, wspace=0.2, hspace=0.75)

    # define x-axis
    t = linspace(-0.5, 2, int(fs*2.5))

    # first column
    ax1 = fig.add_subplot(gs[1:3, 0])
    ax2 = fig.add_subplot(gs[3:5, 0])
    ax3 = fig.add_subplot(gs[5:7, 0])
    ax4 = fig.add_subplot(gs[7:9, 0])
    ax5 = fig.add_subplot(gs[9:11, 0])
    # second column
    ax6 = fig.add_subplot(gs[0:2, 1])
    ax7 = fig.add_subplot(gs[2:4, 1])
    ax8 = fig.add_subplot(gs[4:6, 1])
    ax9 = fig.add_subplot(gs[6:8, 1])
    ax10 = fig.add_subplot(gs[8:10, 1])
    ax11 = fig.add_subplot(gs[10:, 1])
    # third column
    ax12 = fig.add_subplot(gs[1:3, 2])
    ax13 = fig.add_subplot(gs[3:5, 2])
    ax14 = fig.add_subplot(gs[5:7, 2])
    ax15 = fig.add_subplot(gs[7:9, 2])
    ax16 = fig.add_subplot(gs[9:11, 2])

    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10,
            ax11, ax12, ax13, ax14, ax15, ax16]
    titles = ['Sp2', 'Sp5', 'Sp8', 'Sp11', 'Sp14', 'Sp1', 'Sp4', 'Sp7',
              'Sp10', 'Sp13', 'Sp16', 'Sp3', 'Sp6', 'Sp9', 'Sp12', 'Sp15',]
    for idx, ax in enumerate(axes):
            stand_dev_upper = data[idx,:] \
                            + 0.35*std(data_runs[:, idx, :], axis = 0)
            stand_dev_lower = data[idx,:] \
                            - 0.35*std(data_runs[:, idx, :], axis = 0)
            for idx_run, _ in enumerate(data_runs):
                ax.plot(t, data_runs[idx_run, idx, :], 'grey',
                        alpha = 0.25, linewidth = '0.8')
                ax.fill_between(t, stand_dev_lower, stand_dev_upper,
                                alpha = 0.15, color = 'tomato')
            max_val = max(data_runs[:, idx, :])
            max_val = 3
            min_val = min(data_runs[:, idx, :])
            min_val = -3
            ax.vlines(0, min_val, max_val, colors = 'grey',
                      linestyles = 'dotted')
            ax.vlines(1, min_val, max_val, colors = 'grey',
                      linestyles = 'dotted')
            ax.plot(t, data[idx,:], 'k', alpha = 1, linewidth = '0.8')
            ax.set_title(titles[idx], fontweight = 'bold', fontsize = 15)
            
            ax.spines[['right', 'top']].set_visible(False)
            ax.tick_params(axis = 'both', which = 'both', length = 0)
            ax.set_yticks([-2, 0, 2],['-2', '0', '2'], fontsize = 12)
            ax.set_xticks([0, 1, 2],['0', '1', '2'], fontsize = 12)

            ax.set_ylim([-3, 3])
            #plt.yticks([-2.5, 0, 2.5], ['-2.5', '0', '2.5'])

            if not (ax == ax5 or ax == ax11 or ax == ax16):
                 ax.set_xticks([], [])

def ApplyFilters(filter_coeffs, signal_og, fs, t_stim, idx,
                 show_filter_impact = False):
    """
    For a given list fo filter coefficients, this functions applies them
    all to a given signal.

    Input:
        - filter_coeffs: List of tuples that contain the filter coefficients.
        [list, tuple]
        - signal_og: Contains the original signal before filtering.
        [Numpy array]
        - fs: Sampling frequency. [int]
        - t_stim: Stimulation time during the experiment. [float] 
        - idx: indices that refer to that part of the signal that should be
        visualized. [list, int]
        - show_filter_impact: boolean variable to decide wether the filter
        impact should be plotted. [bool]
    
    Output:
        - signal_filtered: Signal after it was filtered with all the given
        filter. [Numpy array]
    """
    idx_run = idx[0]
    idx_chan = idx[1]
    # set filter values
    low_pass = filter_coeffs[0]
    high_pass = filter_coeffs[1]
    notch_filter_35 = filter_coeffs[2]
    notch_filter_50 = filter_coeffs[3]
    
    # filter signals
    signal_filt_1 = filtfilt(low_pass[0], low_pass[1], signal_og)
    signal_filt_2 = filtfilt(high_pass[0], high_pass[1], signal_filt_1)
    signal_filt_3 = filtfilt(notch_filter_35[0],
                                        notch_filter_35[1], signal_filt_2)
    signal_filt_4 = filtfilt(notch_filter_50[0],
                                        notch_filter_50[1],
                                        signal_filt_3)

    # visualize filter impact
    if show_filter_impact:
        if idx_chan == 8 and idx_run == 3:
            rcParams['lines.linewidth'] = 0.8
            rcParams['figure.figsize'] = (12, 8)
            #rcParams['marker.markersize'] = 1
            _, axes = plt.subplots(6,2)

            ax11 = axes[0][0]
            ax12 = axes[0][1]

            ax21 = axes[1][0]
            ax22 = axes[1][1]

            ax31 = axes[2][0]
            ax32 = axes[2][1]

            ax41 = axes[3][0]
            ax42 = axes[3][1]

            ax51 = axes[4][0]
            ax52 = axes[4][1]

            # frequency domain of original signal
            sig_fft = fft.fft(signal_og[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            sig_fft = fft.fftshift(abs(sig_fft))
            # x axis for frequency plots
            freq = linspace(-fs/2, fs/2, len(sig_fft)+1)
            freq = freq[0:-1]
            # x axis for time plots
            t = linspace(-0.5,2.0,int(fs*2.5))

            # show original signal in frequency and time domain
            ax11.plot(t, signal_og[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            ax12.stem(freq,sig_fft, markerfmt = '.')
            ax12.set_xlim([0,150])
            # show signal in frequency and time domain after lowpass
            ax21.plot(t, signal_filt_1[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            sig_fft_1 = signal_filt_1[t_stim-int(fs*0.5):t_stim+int(fs*2)]
            sig_fft_1 = fft.fftshift(abs(fft.fft(sig_fft_1)))
            ax22.stem(freq,sig_fft_1, markerfmt = '.')
            ax22.set_xlim([0,150])
            # show signal in frequency and time domain after
            # lowpass and highpass
            ax31.plot(t, signal_filt_2[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            sig_fft_2 = signal_filt_2[t_stim-int(fs*0.5):t_stim+int(fs*2)]
            sig_fft_2 = fft.fftshift(abs(fft.fft(sig_fft_2)))
            ax32.stem(freq,sig_fft_2, markerfmt = '.')
            ax32.set_xlim([0,150])
            # show signal in frequency and time domain after lowpass,
            # highpass and 35 Hz notch filter
            ax41.plot(t, signal_filt_3[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            sig_fft_3 = signal_filt_3[t_stim-int(fs*0.5):t_stim+int(fs*2)]
            sig_fft_3 = fft.fftshift(abs(fft.fft(sig_fft_3)))
            ax42.stem(freq,sig_fft_3, markerfmt = '.')
            ax42.set_xlim([0,150])
            # show signal in frequency and time domain after lowpass,
            # highpass and 50Hz notch filter
            ax51.plot(t, signal_filt_4[t_stim-int(fs*0.5):t_stim+int(fs*2)])
            sig_fft_4 = signal_filt_4[t_stim-int(fs*0.5):t_stim+int(fs*2)]
            sig_fft_4 = fft.fftshift(abs(fft.fft(sig_fft_4)))
            
            ax52.stem(freq,sig_fft_4, markerfmt = '.')
            ax52.set_xlim([0,150])

            plt.show()

    return signal_filt_4

def FilterData(data, filter_params, marker_exp):
    """
    Functions that creates the necessary filter coefficients
    and applies them by calling the function ApplyFilters().

    Input:
        - data: Data that has to be filtered. [Numpy array]
        - filter_params: Filter parameters, that contain the
        frequencies that for low-, high-pass and notch filter.
        [list, int]
        - marker_exp: Marker data, needed for outlier rejection.
        [Numpy array]
    
    Output:
        data_clean: Processed data. [Numpy array]
    """
    low_pass = filter_params[0]
    high_pass = filter_params[1]
    notch_filter_35 = filter_params[2][0]
    notch_filter_50 = filter_params[2][1]
    fs = filter_params[3]
    show = filter_params[4]
    # create filters
    b1, a1 = butter(N=5, Wn = low_pass, btype = 'highpass', fs = fs,)
    b2, a2 = butter(N=3, Wn = high_pass, btype = 'lowpass', fs = fs)
    b3, a3 = iirnotch(w0 = notch_filter_35,fs = fs, Q = 5)
    b4, a4 = iirnotch(w0 = notch_filter_50,fs = fs, Q = 5)

    filter_coeffs = [[b1, a1], [b2, a2], [b3, a3], [b4, a4]]

    new_order = array([2,5,8,11,14,1,4,7,10,13,16,3,6,9,12,15]) - 1
    data_order = []

    for idx_run, run in enumerate(data):
        new_run = ones_like(run)
        for idx_ch, channel in enumerate(run):
            # apply filters
            new_run[idx_ch] = ApplyFilters(filter_coeffs, channel, fs,
                                           marker_exp[idx_run][15][1],
                                           (idx_run,idx_ch),
                                           show_filter_impact = show)
        data_order.append(new_run[new_order,:])

    # remove peaks / outliers
    # OutlierRejection() functions is missing
    #data_clean = OutlierRejection(data_order, fs, marker_exp)
    print('Data processed.')
    return data_order   # change to data_clean

def AverageData(data, marker, fs):
    """
    With the data and the marker timing this function averages the data.

    Input:
        - data: SPC data of one run of the experiment. [Numpy array]
        - marker: Marker data that was sent during the
        experiment. [Numpy array]
        - fs: Sampling frequency during the experiment. [int]
    
    Output:
        - avg_run: Averaged data of one run. [Numpy array]
    """
    num_runs = len(data)
    runs = zeros([num_runs, 16, int(fs*2.5)])
    for run in range(num_runs):
        
        marker_i = marker[run]
        start_stims = where(marker_i == 3)[0]
        data_i = data[run]

        trials = zeros([start_stims.shape[0], 16, int(fs * 2.5)])

        for trial in range(start_stims.shape[0]):
            trials[trial, :, :] = data_i[:,start_stims[trial] \
                                         - int(fs/2):start_stims[trial] \
                                         + fs * 2]
        runs[run,:,:] = mean(trials, 0)

    avg_run = mean(runs, 0)
    return avg_run

def GetMarker(marker):
    """
    Extracts the timing of the marker data for specific points of the
    experiment, e.g., start_trial, start_stim, ... .

    Input:
        - marker: Marker data that was sent during the
        experiment. [Numpy array]

    Output:
        - marker_exp: Timings of all the times that are of interest.
        [list]
    
    """
    num_runs = len(marker)
    marker_exp = []
    for run in range(num_runs):
        marker_runs = []
        
        marker_i = marker[run]
        
        start_trials = where(marker_i == 2)[0]
        start_stims = where(marker_i == 3)[0]
        start_post_stims = where(marker_i == 4)[0]
        end_trials = where(marker_i == 5)[0]

        if not start_trials.shape == end_trials.shape \
               == start_stims.shape == start_post_stims.shape:
            raise IndexError('Not enough markers were found.')

        for trial in range(start_trials.shape[0]):
            marker_per_trial = [start_trials[trial], start_stims[trial],
                                start_post_stims[trial], end_trials[trial]]
            marker_runs.append(marker_per_trial)
        marker_exp.append(marker_runs)
    print('All markers found.')
    return marker_exp

def ProcessData(data, marker, num_participants, fs, main_path, user_ID):
    """
    This functions is called in the main file. Gets the whole raw data
    structure and returns the filtered and averaged data for each participant
    as well as the overall averaged data.

    Input:
        - data: SPC data of one run of the experiment. [Numpy array]
        - marker: Marker data that was sent during the
        - num_participants: Number of participants in the experiment. [int]
        - fs: Sampling frequency during the experiment. [int]
         - main_path: Contains the path of the main file [pathlib object]
        - user_ID: Is a short ID that points to the folder structure of
        interest. [str]

    Output:
        - data_filt_avg: Processed data that is averages over all trials,
        runs and participants. [Numpy array]
        - processed_data_participants: Processed data that is averages over
        all trials and runs. [Numpy array]
    """
    # define filter params  
    low = 0.5         # low pass filter
    high = 60       # high pass filter
    notch = [35, 50]      # notch filter
    show = False    # show processed data
    filter_params = [low, high, notch, fs, show]

    processed_data_participants = zeros(shape = (num_participants,16,1280))
    for i in range(num_participants):
        # get marker
        marker_exp = GetMarker(marker[i])
        # filter data
        data_filt = FilterData(data[i], filter_params, marker_exp)
        # avergae data
        data_filt_avg = AverageData(data_filt, marker[i], fs)
        processed_data_participants[i,:,:] = data_filt_avg
    
    data_filt_avg = mean(processed_data_participants, 0)
    #Create16ChannelPlot(data_filt_avg, processed_data_participants, fs)

    # save processed data
    user_folder = 'user_' + user_ID
    data_avg_path = main_path / 'data' / user_folder / 'processed_data' \
                    / 'data_avg.npy'
    data_runs_path = main_path / 'data' / user_folder /'processed_data' \
                    / 'data_participants.npy'
    save(data_avg_path, data_filt_avg)
    save(data_runs_path, data_filt_avg)
    
    
    return data_filt_avg, processed_data_participants

if __name__ == '__main__':
    ...
