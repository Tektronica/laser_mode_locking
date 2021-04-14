import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy.signal import hilbert

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

    if windfunc == 'rectangular':
        M = int(fs / ldf)
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


def get_random_phase(dpi=np.pi / 2):
    return np.random.randint(5) * dpi


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
    return yt_list


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
    print(lowermin, uppermin)
    return lowermin, uppermin


def get_envelope_FWHM(envelope, fs):
    """
    computes the Full-Wave, Half-Maximum of the envelope
    """
    # skips the first peak. Data set potentially starts at the first peak, which prevents left min to be determined.
    length = envelope.size
    peak = int(np.argmax(envelope[int(length / 4):])+int(length / 4))
    print(length)
    print(int(np.argmax(envelope[int(length / 4):])))
    print(int(length / 4))
    print('peak index', peak)
    print('envelope', envelope[peak])
    lowermin, uppermin = find_range(envelope, peak)

    # fwhm_index = np.where(np.isclose(envelope[lowermin:uppermin], envelope[peak]/2, atol=980e-6))
    # fwhm_val = envelope[fwhm_index[0]]
    # fwhm_width = np.diff(fwhm_index)/fs

    # fwhm_index = np.where(np.isclose(envelope[lowermin:uppermin], envelope[peak]/2, atol=980e-6))
    band = uppermin - lowermin
    fwhm_band = int(band/2)
    fwhm_val = envelope[lowermin + int(fwhm_band/2)]
    print('band',band)
    print('fwhm', fwhm_val)
    fwhm_width = fwhm_band/fs

    return fwhm_val, fwhm_width


def rms_flat(a):
    """
    Return the root mean square of all the elements of *a*, flattened out.
    """
    return np.sqrt(np.mean(np.absolute(a) ** 2))


def simulation():
    # http://www.uobabylon.edu.iq/eprints/publication_2_14877_1775.pdf

    # PARAMETERS -------------------------------------------------------------------------------------------------------
    WINDOW_FUNC = 'rectangular'
    fc = 473.613e12  # actual vacuum frequency of HeNe (632.991 nm)
    error = 0.01
    emitted_modes = 7  # number of modes/tones
    n = 1.0  # index of refraction
    # laser_bw = 1.5e9  # HeNe
    laser_bw = fc * 0.1
    BANDWIDTH_SHAPE = 'flat-top'  # gaussian
    print((fc-laser_bw/2)/1e12, (fc+laser_bw/2)/1e12)
    random_phase = False
    gaussian_profile = 20  # Gaussian profile standard deviation

    print('laser bandwidth:', laser_bw / 1e12, 'THz')

    FWHM = laser_bw
    print('full wave, half maximum:', laser_bw / 1e12, 'THz')

    wavelength = SPEED_OF_LIGHT / fc
    print('wavelength, lambda:', round(wavelength * 1e9, 3), 'nm')
    print()

    df_max = laser_bw / emitted_modes
    print('max frequency separation for number of emitted modes, df:', round(df_max / 1e9, 3), 'GHz')

    cavity_modes = np.ceil(fc / (n * df_max))
    print('cavity modes, m:', cavity_modes)

    cavity_length = cavity_modes * wavelength / 2
    print('cavity length, L:', round(cavity_length * 1e2, 3), 'cm', round(cavity_length * 1e3, 3), '(mm)')

    cavity_df = SPEED_OF_LIGHT / (2 * n * cavity_length)
    print('frequency separation of cavity, df:', round(cavity_df / 1e6, 3), 'MHz')

    longitudinal_modes = int(laser_bw / cavity_df)  # the number of modes supported by the laser bandwidth
    print('longitudinal modes supported:', longitudinal_modes)
    print()

    laser_power = 0.0
    print('laser power:', round(laser_power / 1e3, 3), 'mW')

    # TIME BASE --------------------------------------------------------------------------------------------------------
    fs = fc * 100
    main_lobe_error = min(cavity_df / (50 * fc), error)

    N = getWindowLength(f0=fc, fs=fs, windfunc=WINDOW_FUNC, error=main_lobe_error)

    runtime = N / fs
    print('runtime:', round(runtime * 1e12, 3), 'ps')
    N_range = np.arange(0, N, 1)
    xt = N_range / fs

    # WAVEFORM GENERATOR -----------------------------------------------------------------------------------------------
    yt_list = get_waveforms(xt, fc,
                            df=cavity_df, modes=emitted_modes,
                            bw=laser_bw, bw_shape=BANDWIDTH_SHAPE, random_phase=random_phase)

    yt = np.sum(yt_list, axis=0)  # element-wise summation
    yt_envelope = compute_envelope(yt)
    fwhm_val, fwhm_width = get_envelope_FWHM(yt_envelope, fs)
    print()
    print('FWHM value:', round(fwhm_val, 2))
    print('FWHM width:', round(fwhm_width*1e12, 3), 'ps')

    xf_fft, yf_fft, xf_rfft, yf_rfft, fft_length, main_lobe_width = windowed_fft(yt, fs, N, WINDOW_FUNC)
    yf_smooth = get_gaussian(xf_rfft, fc, laser_bw, bw_shape=BANDWIDTH_SHAPE)

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

    dim_left = (fc - laser_bw / 2) / xf_scale
    dim_right = (fc + laser_bw / 2) / xf_scale
    dim_height = get_gaussian(dim_left * xf_scale, fc, laser_bw)

    # Arrow dimension line update ----------------------------------------------------------------------------------
    # https://stackoverflow.com/a/48684902 -------------------------------------------------------------------------
    arrow_dim_obj.xy = (dim_left, dim_height)
    arrow_dim_obj.set_position((dim_right, dim_height))
    arrow_dim_obj.textcoords = ax2.transData

    # dimension text update ----------------------------------------------------------------------------------------
    dim_text.set_position((dim_left + ((laser_bw / xf_scale) / 2), dim_height))
    dim_text.set_text(f"FWHM: {round(FWHM / 1e12, 3)} THz")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    simulation()
