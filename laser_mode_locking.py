import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy.signal import hilbert
import math

pylab_params = {'legend.fontsize': 'medium',
                'font.family': 'Segoe UI',
                'axes.titleweight': 'bold',
                'figure.figsize': (15, 5),
                'axes.labelsize': 'medium',
                'axes.titlesize': 'medium',
                'xtick.labelsize': 'medium',
                'ytick.labelsize': 'medium'}
pylab.rcParams.update(pylab_params)

SPEED_OF_LIGHT = 299792458  # speed of light (3e8 m/s)


def getWindowLength(f0=10e3, fs=2.5e6, windfunc='blackman', mlw=0.01, mainlobe_type='relative'):
    """
    Computes the window length of the measurement. An error is expressed since the main lobe width is directly
    proportional to the number of cycles captured. The minimum value of M correlates to the lowest detectable frequency
    by the windowing function. For instance, blackman requires a minimum of 6 period cycles of the frequency of interest
    in order to express content of that lobe in the DFT. Sampling frequency does not play a role in the width of the
    lobe, only the resolution of the lobe.

    :param mlw: Mainlobe width (can be specified relative or absolute w.r.t. fundamental)
    :param mainlobe_type: Mainlobe width can be set relative to the signal frequency or as an absolute width
    independent of signal frequency

    :param f0: fundamental frequency of signal
    :param fs: sampling frequency
    :param windfunc: "Rectangular", "Bartlett", "Hanning", "Hamming", "Blackman"
    :return: window length of integer value (number of time series
    samples collected)
    """
    # lowest detectable frequency by window
    # aka - the main lobe width
    if mainlobe_type == 'relative':
        ldf = f0 * mlw
    elif mainlobe_type == 'absolute':
        ldf = mlw
    else:
        raise ValueError('Incorrect main lobe type used!\nSelection should either be relative or absolute.')

    if windfunc == 'rectangular':
        M = int(1 * (fs / ldf))
    elif windfunc in ('bartlett', 'hanning', 'hamming'):
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
    if windfunc == 'rectangular':
        w = np.ones(N)
        main_lobe_width = 2 * (Fs / N)
    elif windfunc == 'bartlett':
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


def get_random_phase(dpi=np.pi / 6, distribution='gaussian'):
    if distribution == 'random':
        return np.random.randint(5) * dpi
    else:
        return int(np.random.normal(5)) * dpi


def get_gaussian(f, fc, bw, bw_shape='gaussian'):
    if bw_shape == 'gaussian':
        return np.exp(-0.5 * ((f - fc) / (0.849322 * bw / 2)) ** 2)
    else:
        return np.where(f > (fc - bw / 2), np.where(f <= (fc + bw / 2), 1, 0), 0)


def get_waveforms(xt, fc, df, modes=5, bw=1e9, bw_shape='gaussian', random_phase=False):
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
        A = get_gaussian(f, fc, bw, bw_shape=bw_shape)
        yt_list[mode] = A * np.sin(2 * np.pi * f * xt + phase)

    return np.array(yt_list)


def compute_envelope(signal):
    # returns absolute value of the hilbert transformation
    return abs(hilbert(signal))


def find_range(f, x):
    """
    Find range between nearest local minima from peak at index x
    """
    lowermin = 0.0
    uppermin = 0.0

    for i in np.arange(x + 1, len(f)):
        if f[i + 1] >= f[i]:
            uppermin = i
            break

    for i in np.arange(x - 1, 0, -1):
        if f[i] <= f[i - 1]:
            lowermin = i + 1
            break

    return lowermin, uppermin


def get_envelope_FWHM(envelope, fs):
    """
    computes the Full-Wave, Half-Maximum of the envelope
    """
    # skips the first peak. Data set potentially starts at the first peak, which prevents left min to be determined.
    length = envelope.size
    chunk = int(length / 4)
    peak_index = int(np.argmax(envelope[chunk:3*chunk]) + chunk)
    lowermin, uppermin = find_range(envelope, peak_index)
    fwhm_val = envelope[peak_index]/2

    pulse = envelope[lowermin:uppermin]
    fwhm_index = np.where(np.isclose(pulse, fwhm_val, atol=fwhm_val/(pulse.size/2)))
    fwhm_width = np.diff(fwhm_index)/fs

    try:
        return fwhm_val, fwhm_width[0][0]
    except IndexError:
        print('index was out of bounds. size 0 more than likely.')
        return fwhm_val, np.NaN


def rms_flat(a):
    """
    Return the root mean square of all the elements of *a*, flattened out.
    """
    # https://stackoverflow.com/a/17463210
    # https://code.activestate.com/recipes/393090/
    # https://stackoverflow.com/a/33004170
    sqr = np.absolute(a) ** 2
    mean = math.fsum(sqr) / len(sqr)  # computed from partial sums
    return np.sqrt(mean)


def worker(params):
    fc, laser_bw, emitted_modes, refraction_index, bandwidth_shape, window, MLW, random_phase = params

    gaussian_profile = 20  # Gaussian profile standard deviation
    fc = fc * 1e12
    laser_bw = fc * laser_bw

    wavelength = SPEED_OF_LIGHT / fc
    df_max = laser_bw / emitted_modes
    cavity_modes = np.ceil(fc / (refraction_index * df_max))
    cavity_length = cavity_modes * wavelength / 2
    cavity_df = SPEED_OF_LIGHT / (2 * refraction_index * cavity_length)
    longitudinal_modes = int(laser_bw / cavity_df)  # the number of modes supported by the laser bandwidth

    # TIME BASE --------------------------------------------------------------------------------------------------------
    fs = fc * 100
    main_lobe_error = min(cavity_df / (50 * fc), MLW)

    N = getWindowLength(f0=fc, fs=fs, windfunc=window, mlw=main_lobe_error, mainlobe_type='relative')

    runtime = N / fs

    N_range = np.arange(0, N, 1)
    xt = N_range / fs

    # WAVEFORM GENERATOR -----------------------------------------------------------------------------------------------
    yt_list = get_waveforms(xt, fc,
                            df=cavity_df, modes=emitted_modes,
                            bw=laser_bw, bw_shape=bandwidth_shape, random_phase=random_phase)

    yt = np.sum(yt_list, axis=0)  # element-wise summation
    yt_envelope = compute_envelope(yt)
    fwhm_val, fwhm_width = get_envelope_FWHM(yt_envelope, fs)

    xf_fft, yf_fft, xf_rfft, yf_rfft, fft_length, main_lobe_width = windowed_fft(yt, fs, N, window)
    yf_smooth = get_gaussian(xf_rfft, fc, laser_bw, bw_shape=bandwidth_shape)

    # LIMITS -----------------------------------------------------------------------------------------------------------
    xt1_left = 0  # show the start of the data
    xt1_right = 2 / cavity_df  # want to display two periods of the frequency step

    xt2_left = 0
    xt2_right = 4 / cavity_df  # what to display two periods of the frequency step

    bound = 2 * cavity_df * (emitted_modes - 1) / 2
    xf_left = max(0.0, fc - bound)
    xf_right = min(fc + bound, fs / 2)

    # ANNOTATIONS  -----------------------------------------------------------------------------------------------------
    xf_scale = 1e12
    dim_left = (fc - laser_bw / 2) / xf_scale
    dim_right = (fc + laser_bw / 2) / xf_scale
    dim_height = get_gaussian(dim_left * xf_scale, fc, laser_bw)

    dim_label = f"FWHM: {round(laser_bw / 1e12, 3)} THz"
    dim_label_pos = (laser_bw / xf_scale) / 2

    data = (wavelength, laser_bw, df_max, cavity_modes, cavity_length, cavity_df, longitudinal_modes, fwhm_val, fwhm_width, runtime)
    plot_data = (xt, yt_list, yt, yt_envelope, xf_rfft, yf_rfft, yf_smooth)
    plot_limits = xt1_left, xt1_right, xt2_left, xt2_right, xf_left, xf_right, dim_left, dim_right, dim_height, dim_label, dim_label_pos

    return data, plot_data, plot_limits


def simulation():
    # http://www.uobabylon.edu.iq/eprints/publication_2_14877_1775.pdf

    # PARAMETERS -------------------------------------------------------------------------------------------------------
    window = 'rectangular'
    fc = 473.613  # (THz) actual vacuum frequency of HeNe (632.991 nm)
    MLW = 0.01  # relative
    emitted_modes = 15  # number of modes/tones
    refraction_index = 1.0  # index of refraction
    # laser_bw = 1.5e9  # HeNe
    laser_bw = 0.1
    bandwidth_shape = 'flat-top'  # gaussian
    random_phase = True

    params = (fc, laser_bw, emitted_modes, refraction_index, bandwidth_shape, window, MLW, random_phase)
    data, plot_data, plot_limits = worker(params)

    wavelength, laser_bw, df_max, cavity_modes, cavity_length, cavity_df, longitudinal_modes, fwhm_val, fwhm_width, runtime = data
    xt1_left, xt1_right, xt2_left, xt2_right, xf_left, xf_right, dim_left, dim_right, dim_height, dim_label, dim_label_pos = plot_limits

    fc = fc * 1e12
    FWHM = laser_bw
    fs = fc * 100

    print('laser bandwidth:', laser_bw / 1e12, 'THz')
    print('full wave, half maximum:', laser_bw / 1e12, 'THz')
    print('wavelength, lambda:', round(wavelength * 1e9, 3), 'nm')
    print()
    print('max frequency separation for number of emitted modes, df:', round(df_max / 1e9, 3), 'GHz')
    print('cavity modes, m:', cavity_modes)
    print('cavity length, L:', round(cavity_length * 1e2, 3), 'cm', round(cavity_length * 1e3, 3), '(mm)')
    print('frequency separation of cavity, df:', round(cavity_df / 1e9, 3), 'GHz')
    print('longitudinal modes supported:', longitudinal_modes)
    print()

    xt, yt_list, yt, yt_envelope, xf_rfft, yf_rfft, yf_smooth = plot_data

    # PLOT GENERATOR -----------------------------------------------------------------------------------------------
    figure = plt.figure(figsize=(12.8, 9.6), constrained_layout=False)  # default: figsize=(6.4, 4.8)
    ax1 = figure.add_subplot(311)
    ax2 = figure.add_subplot(312)
    ax3 = figure.add_subplot(313)

    xt_scale = 1e12
    xf_scale = 1e12

    for yt_data in yt_list:
        temporal1, = ax1.plot(xt * xt_scale, yt_data, '-')  # All signals
    temporal2, = ax2.plot(xt * xt_scale, yt, '-')  # The summation of all signals
    temporal3, = ax2.plot(xt * xt_scale, yt_envelope, '-')  # The envelope of the summation
    spectral1, = ax3.plot(xf_rfft / xf_scale, np.abs(yf_rfft), '-', color='#C02942')  # The spectral plot of sum
    spectral2, = ax3.plot(xf_rfft / xf_scale, yf_smooth, '-')  # The spectral plot of sum

    # LIMITS -----------------------------------------------------------------------------------------------------------
    xt1_left = 0  # show the start of the data
    xt1_right = 2 / cavity_df  # want to display two periods of the frequency step
    ax1.set_xlim(left=xt1_left * xt_scale, right=xt1_right * xt_scale)
    # ax1.set_ylim(bottom=yf_btm, top=yf_top)

    xt2_left = 0
    xt2_right = 6 / cavity_df  # what to display two periods of the frequency step
    ax2.set_xlim(left=xt2_left * xt_scale, right=xt2_right * xt_scale)
    # ax2.set_ylim(bottom=yf_btm, top=yf_top)

    bound = 2 * cavity_df * (emitted_modes - 1) / 2
    xf_left = max(0.0, fc - bound)
    xf_right = min(fc + bound, fs / 2)
    ax3.set_xlim(left=xf_left / xf_scale, right=xf_right / xf_scale)
    # ax3.set_ylim(bottom=yf_btm, top=yf_top)

    ax1.set_title('SAMPLED TIMED SERIES DATA')
    ax1.set_xlabel('TIME (ps)')
    ax1.set_ylabel('AMPLITUDE')

    ax2.set_title('SUMMATION OF ALL MODES/TONES')
    ax2.set_xlabel('TIME (ps)')
    ax2.set_ylabel('AMPLITUDE')

    ax3.set_title('SPECTRAL DATA')
    ax3.set_xlabel('FREQUENCY (THz)')
    ax3.set_ylabel('MAGNITUDE (V)')

    # Annotations --------------------------------------------------------------------------------------------------
    arrow_dim_obj = ax3.annotate("", xy=(0, 0), xytext=(0, 0),
                                 textcoords=ax3.transData, arrowprops=dict(arrowstyle='<->'))
    bbox = dict(fc="white", ec="none")
    dim_text = ax3.text(0, 0, "", ha="center", va="center", bbox=bbox)

    # Arrow dimension line update ----------------------------------------------------------------------------------
    # https://stackoverflow.com/a/48684902 -------------------------------------------------------------------------
    arrow_dim_obj.xy = (dim_left, dim_height)
    arrow_dim_obj.set_position((dim_right, dim_height))
    arrow_dim_obj.textcoords = ax2.transData

    # dimension text update ----------------------------------------------------------------------------------------
    dim_text.set_position((dim_left + dim_label_pos, dim_height))
    dim_text.set_text(dim_label)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    simulation()
