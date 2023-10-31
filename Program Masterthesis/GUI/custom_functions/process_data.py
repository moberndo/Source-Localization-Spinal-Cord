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
    #print('All markers found.')
    return marker_exp

def ProcessData(data, marker, num_participants, fs, main_path):
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
    
    
    return data_filt_avg, processed_data_participants

if __name__ == '__main__':
    ...
