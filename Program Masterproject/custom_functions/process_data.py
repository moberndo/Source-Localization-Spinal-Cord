""" ####################################################################### 
    Title: Construction of an FE model of the spinal cord to simulate the voltage distribution caused by a dipole
    Subtitle: Process data Function
    Author: Markus E. Oberndorfer
####################################################################### """

# imports
from scipy import fft
from scipy.signal import filtfilt, butter, iirnotch
from numpy import ones_like, array, argsort, abs, linspace
import matplotlib.pyplot as plt

def ApplyFilters(filter_coeffs, signal):
    for filter_coeff in filter_coeffs:
        signal = filtfilt(filter_coeff[0], filter_coeff[1], signal)

    return signal


def ProcessData(data, low = 5, high = 60, notch = 35, fs = 512, show = False):
    # create filters
    b1, a1 = butter(N=5, Wn = low, btype = 'highpass', fs = fs,)
    b2, a2 = butter(N=3, Wn = high, btype = 'lowpass', fs = fs)
    b3, a3 = iirnotch(w0 = notch,fs = fs, Q = 30)

    filter_coeffs = [[b1, a1], [b2, a2], [b3, a3]]


    new_order = array([2,5,8,11,14,1,4,7,10,13,16,3,6,9,12,15]) - 1
    new_order_ix = argsort(new_order)
    #data_order = ones_like(data)
    data_order = []

    for idx_r, run in enumerate(data):
        new_run = ones_like(run)
        for idx_ch, channel in enumerate(run):
            # apply filters
            new_run[idx_ch] = ApplyFilters(filter_coeffs, channel)
        data_order.append(new_run[new_order,:])

    print('Data processed.')
    if show == True:
        sig_fft = fft.fft(data_order[0][0,:])
        freq = linspace(-fs/2, fs/2, len(sig_fft)+1)
        freq = freq[0:-1]

        sig_fft2 = fft.fft(data_order[0][0,:])

        _, ax = plt.subplots()
        ax.stem(freq,fft.fftshift(abs(sig_fft)))

        _, ax2 = plt.subplots()
        ax2.stem(freq,fft.fftshift(abs(sig_fft2)))

    
    return data_order

