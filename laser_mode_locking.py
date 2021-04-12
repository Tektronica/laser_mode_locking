import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

pylab_params = {'legend.fontsize': 'medium',
                'font.family': 'Segoe UI',
                'axes.titleweight': 'bold',
                'figure.figsize': (15, 5),
                'axes.labelsize': 'medium',
                'axes.titlesize': 'medium',
                'xtick.labelsize': 'medium',
                'ytick.labelsize': 'medium'}
pylab.rcParams.update(pylab_params)


def getWindowLength(f0=10e3, fs=2.5e6, windfunc='blackman', error=0.1):
    """
    Computes the window length of the measurement. An error is expressed since the main lobe width is directly
    proportional to the number of cycles captured. The minimum value of M correlates to the lowest detectable frequency
    by the windowing function. For instance, blackman requires a minimum of 6 period cycles of the frequency of interest
    in order to express content of that lobe in the DFT. Sampling frequency does not play a role in the width of the
    lobe, only the resolution of the lobe.

    :param f0: frequency of interest
    :param fs: sampling frequency
    :param windfunc: "Rectangular", "Bartlett", "Hanning", "Hamming", "Blackman"
    :param error: 100% error suggests the lowest detectable frequency is the fundamental
    :return: window length of integer value (number of time series samples collected)
    """
    # lowest detectable frequency by window
    ldf = f0 * error

    if windfunc == 'Rectangular':
        M = int(fs / ldf)
    elif windfunc in ('Bartlett', 'Hanning', 'Hamming'):
        M = int(4 * (fs / ldf))
    elif windfunc == 'blackman':
        M = int(6 * (fs / ldf))
    else:
        raise ValueError('Not a valid windowing function.')

    return M


def windowed_fft(yt, Fs, N, windfunc='blackman'):
    # remove DC offset
    yt -= np.mean(yt)

    # Calculate windowing function and its length ----------------------------------------------------------------------
    if windfunc == 'bartlett':
        w = np.bartlett(N)
        main_lobe_width = 4 * (Fs / N)
    elif windfunc == 'hanning':
        w = np.hanning(N)
        main_lobe_width = 4 * (Fs / N)
    elif windfunc == 'hamming':
        w = np.hamming(N)
        main_lobe_width = 4 * (Fs / N)
    elif windfunc == 'blackman':
        w = np.blackman(N)
        main_lobe_width = 6 * (Fs / N)
    else:
        # TODO - maybe include kaiser as well, but main lobe width varies with alpha
        raise ValueError("Invalid windowing function selected!")

    # Calculate amplitude correction factor after windowing ------------------------------------------------------------
    # https://stackoverflow.com/q/47904399/3382269
    amplitude_correction_factor = 1 / np.mean(w)

    # Calculate the length of the FFT ----------------------------------------------------------------------------------
    if (N % 2) == 0:
        # for even values of N: FFT length is (N / 2) + 1
        fft_length = int(N / 2) + 1
    else:
        # for odd values of N: FFT length is (N + 1) / 2
        fft_length = int((N + 2) / 2)

    """
    Compute the FFT of the signal Divide by the length of the FFT to recover the original amplitude. Note dividing 
    alternatively by N samples of the time-series data splits the power between the positive and negative sides. 
    However, we are only looking at one side of the FFT.
    """
    yf_fft = (np.fft.fft(yt * w) / fft_length) * amplitude_correction_factor

    yf_rfft = yf_fft[:fft_length]
    xf_fft = np.linspace(0.0, Fs, N)
    xf_rfft = np.linspace(0.0, Fs / 2, fft_length)

    return xf_fft, yf_fft, xf_rfft, yf_rfft, fft_length, main_lobe_width


def get_random_phase(dpi=np.pi / 2):
    return np.random.randint(5) * dpi


def get_waveforms(xt, fc, df, modes=7, random_phase=False):
    yt_list = [None] * modes
    for mode in range(modes):

        # calculate location of new spectral component
        new_freq = df * (((mode % 2) * -2) + 1) * int((mode + 1) / 2)
        f = fc + new_freq

        # calculate a random phase if necessary
        phase = 0
        if random_phase:
            # apply a random phase
            phase = get_random_phase()

        # https://www.edmundoptics.com/knowledge-center/application-notes/optics/why-use-a-flat-top-laser-beam/
        # https://micro.magnet.fsu.edu/primer/java/lasers/gainbandwidth/index.html
        # sigma = 1 / (np.sqrt(2 * np.pi))
        # pos = fc/2
        # A = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((pos - 1) / sigma)**2)
        # TODO: calculate bandwidth
        sigma = df*(modes/4)
        A = np.exp(-0.5 * ((f - fc+1) / sigma) ** 2)
        yt_list[mode] = A * np.sin(2 * np.pi * f * xt + phase)
    return yt_list


def rms_flat(a):
    """
    Return the root mean square of all the elements of *a*, flattened out.
    """
    return np.sqrt(np.mean(np.absolute(a) ** 2))


def simulation():
    # http://www.uobabylon.edu.iq/eprints/publication_2_14877_1775.pdf

    # PARAMETERS -------------------------------------------------------------------------------------------------------
    fc = 1000  # center frequency
    df = 10  # spectral separation of modes/tones
    number_of_modes = 21  # number of modes/tones
    gaussian_profile = 20  # Gaussian profile standard deviation
    random_phase = False

    # TIME BASE --------------------------------------------------------------------------------------------------------
    fs = 100000
    main_lobe_error = 2 / 200
    WINDOW_FUNC = 'blackman'
    N = getWindowLength(f0=fc, fs=fs, windfunc=WINDOW_FUNC, error=main_lobe_error)

    runtime = N / fs
    N_range = np.arange(0, N, 1)
    xt = N_range / fs

    # WAVEFORM GENERATOR -----------------------------------------------------------------------------------------------
    yt_list = get_waveforms(xt, fc, df, modes=number_of_modes, random_phase=random_phase)

    # https://mathworld.wolfram.com/FourierTransformGaussian.html
    # https://math.stackexchange.com/questions/1267007/inverse-fourier-transform-of-gaussian
    # https://www.youtube.com/watch?v=a8S26iExR7A
    # sigma = (1/np.sqrt(2*np.pi))
    # gaussian = (1/(sigma * np.sqrt(2*np.pi)))*np.exp((xt/sigma)**2)

    yt = np.sum(yt_list, axis=0)  # element-wise summation
    # ytt = np.convolve(yt, gaussian)

    xf_fft, yf_fft, xf_rfft, yf_rfft, fft_length, main_lobe_width = windowed_fft(yt, fs, N, WINDOW_FUNC)

    # PLOT GENERATOR -----------------------------------------------------------------------------------------------
    figure = plt.figure(figsize=(12.8, 9.6), constrained_layout=True)  # default: figsize=(6.4, 4.8)
    ax1 = figure.add_subplot(311)
    ax2 = figure.add_subplot(312)
    ax3 = figure.add_subplot(313)

    for yt_data in yt_list:
        temporal1, = ax1.plot(xt, yt_data, '-')  # All signals
    temporal2, = ax2.plot(xt, yt, '-')  # The summation of all signals
    spectral1, = ax3.plot(xf_rfft, np.abs(yf_rfft), '-')  # The spectral plot of sum

    # LIMITS -----------------------------------------------------------------------------------------------------------
    xt1_left = 0  # show the start of the data
    xt1_right = 2 / df  # want to display two periods of the frequency step
    ax1.set_xlim(left=xt1_left, right=xt1_right)
    # ax1.set_ylim(bottom=yf_btm, top=yf_top)

    xt2_left = 0
    xt2_right = 6 / df  # what to display two periods of the frequency step
    ax2.set_xlim(left=xt2_left, right=xt2_right)
    # ax2.set_ylim(bottom=yf_btm, top=yf_top)

    bound = 2 * df * (number_of_modes - 1) / 2
    xf_left = max(0.0, fc - bound)
    xf_right = min(fc + bound, fs / 2)
    ax3.set_xlim(left=xf_left, right=xf_right)
    # ax3.set_ylim(bottom=yf_btm, top=yf_top)

    ax1.set_title('SAMPLED TIMED SERIES DATA')
    ax1.set_xlabel('TIME (s)')
    ax1.set_ylabel('AMPLITUDE')

    ax2.set_title('SUMMATION OF ALL MODES/TONES')
    ax2.set_xlabel('TIME (s)')
    ax2.set_ylabel('AMPLITUDE')

    ax3.set_title('SPECTRAL DATA')
    ax3.set_xlabel('FREQUENCY (kHz)')
    ax3.set_ylabel('MAGNITUDE (V)')
    plt.show()


if __name__ == "__main__":
    simulation()
