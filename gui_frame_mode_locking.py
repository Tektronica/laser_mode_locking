import numpy as np

import laser_mode_locking as worker
import wx

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar


# CONFIG ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

# CONFIG \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def to_float(string_val):
    try:
        float_val = float(string_val)
    except ValueError:
        print(f'[ERROR] Invalid Input. Could not convert {string_val} to float.')
        raise ValueError(f'Invalid Input. Could not convert {string_val} to float.')
    else:
        return float_val


def to_integer(string_val):
    try:
        int_val = int(string_val)
    except ValueError:
        print(f'[ERROR] Invalid Input. Could not convert {string_val} to integer.')
        raise ValueError(f'Invalid Input. Could not convert {string_val} to integer.')
    else:
        return int_val


class MyDemoPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, wx.ID_ANY)

        self.frame = parent
        self.left_panel = wx.Panel(self, wx.ID_ANY)
        self.plot_panel = wx.Panel(self, wx.ID_ANY, style=wx.SIMPLE_BORDER)

        # PLOT Panel ---------------------------------------------------------------------------------------------------
        self.figure = plt.figure(figsize=(1, 1))  # look into Figure((5, 4), 75)
        self.canvas = FigureCanvas(self.plot_panel, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()

        self.ax1 = self.figure.add_subplot(311)
        self.ax2 = self.figure.add_subplot(312)
        self.ax3 = self.figure.add_subplot(313)

        self.temporal, = self.ax1.plot([], [], linestyle='-')
        self.temporal_sum, = self.ax2.plot([], [], linestyle='-', marker='')
        self.temporal_hilbert, = self.ax2.plot([], [], linestyle='-', marker='')
        self.spectral, = self.ax3.plot([], [], color='#C02942')
        self.spectral_envelope, = self.ax3.plot([], [], color='tab:blue')

        self.combo_window = wx.ComboBox(self.left_panel, wx.ID_ANY,
                                        choices=["Rectangular", "Bartlett", "Hanning", "Hamming", "Blackman"],
                                        style=wx.CB_DROPDOWN | wx.CB_READONLY)
        self.combo_bandwidth_shape = wx.ComboBox(self.left_panel, wx.ID_ANY,
                                                 choices=["Flat-Top", "Gaussian"],
                                                 style=wx.CB_DROPDOWN | wx.CB_READONLY)

        self.text_ctrl_fc = wx.TextCtrl(self.left_panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER)
        self.text_ctrl_laser_bw = wx.TextCtrl(self.left_panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER)
        self.text_ctrl_mainlobe = wx.TextCtrl(self.left_panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER)
        self.text_ctrl_mainlobe.SetToolTip("Mainlobe width")
        self.text_ctrl_emitted_modes = wx.TextCtrl(self.left_panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER)
        self.text_ctrl_index = wx.TextCtrl(self.left_panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER)
        self.checkbox_random_phase = wx.CheckBox(self.left_panel, wx.ID_ANY, "Random Phase")

        self.report_runtime = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_laser_bw = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_wavelength = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_cavity_modes = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_cavity_length = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_df = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_longitudinal_modes = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_fwhm = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)
        self.report_fwhm_width = wx.TextCtrl(self.left_panel, wx.ID_ANY, "", style=wx.TE_READONLY)

        on_update = lambda event: self.update(event)
        self.Bind(wx.EVT_TEXT_ENTER, on_update, self.text_ctrl_fc)
        self.Bind(wx.EVT_TEXT_ENTER, on_update, self.text_ctrl_laser_bw)
        self.Bind(wx.EVT_TEXT_ENTER, on_update, self.text_ctrl_mainlobe)
        self.Bind(wx.EVT_TEXT_ENTER, on_update, self.text_ctrl_emitted_modes)
        self.Bind(wx.EVT_TEXT_ENTER, on_update, self.text_ctrl_index)
        self.Bind(wx.EVT_COMBOBOX_CLOSEUP, on_update, self.combo_bandwidth_shape)
        self.Bind(wx.EVT_CHECKBOX, on_update, self.checkbox_random_phase)

        self.__set_properties()
        self.__do_layout()
        self.__do_plot_layout()
        self.update(wx.Event)

    def __set_properties(self):
        self.SetBackgroundColour(wx.Colour(240, 240, 240))
        self.canvas.SetMinSize((700, 490))

        self.combo_window.SetSelection(0)
        self.combo_bandwidth_shape.SetSelection(0)
        self.checkbox_random_phase.SetValue(0)

        self.text_ctrl_fc.SetValue("473.613")  # (THz)
        self.text_ctrl_laser_bw.SetValue("0.1")
        self.text_ctrl_index.SetValue("1.0")
        self.text_ctrl_emitted_modes.SetValue("15")
        self.text_ctrl_mainlobe.SetValue("0.01")

        self.report_runtime.SetValue("--")
        self.report_laser_bw.SetValue("--")
        self.report_wavelength.SetValue("--")
        self.report_cavity_modes.SetValue("--")
        self.report_cavity_length.SetValue("--")
        self.report_df.SetValue("--")
        self.report_longitudinal_modes.SetValue("--")
        self.report_fwhm.SetValue("--")
        self.report_fwhm_width.SetValue("--")

    def __do_layout(self):
        sizer_2 = wx.GridSizer(1, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(1, 2, 0, 0)
        grid_sizer_plot = wx.GridBagSizer(0, 0)
        grid_sizer_left_panel = wx.GridBagSizer(0, 0)

        # LEFT PANEL ---------------------------------------------------------------------------------------------------
        # TITLE --------------------------------------------------------------------------------------------------------
        row = 0
        label_1 = wx.StaticText(self.left_panel, wx.ID_ANY, "LASER MODE-LOCKING")
        label_1.SetFont(wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, 0, ""))
        grid_sizer_left_panel.Add(label_1, (row, 0), (1, 3), wx.LEFT | wx.RIGHT | wx.TOP, 5)

        row += 1
        static_line_1 = wx.StaticLine(self.left_panel, wx.ID_ANY)
        static_line_1.SetMinSize((300, 2))
        grid_sizer_left_panel.Add(static_line_1, (row, 0), (1, 3), wx.BOTTOM | wx.RIGHT | wx.TOP, 5)

        # PARAMETERS ---------------------------------------------------------------------------------------------------
        row += 2
        lbl_settings = wx.StaticText(self.left_panel, wx.ID_ANY, "Parameters")
        lbl_settings.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, 0, ""))
        grid_sizer_left_panel.Add(lbl_settings, (row, 0), (1, 3), wx.LEFT | wx.RIGHT, 5)

        row += 1
        static_line_2 = wx.StaticLine(self.left_panel, wx.ID_ANY)
        static_line_2.SetMinSize((300, 2))
        grid_sizer_left_panel.Add(static_line_2, (row, 0), (1, 3), wx.BOTTOM | wx.RIGHT | wx.TOP, 5)

        row += 1
        lbl_fc = wx.StaticText(self.left_panel, wx.ID_ANY, "Fc:")
        grid_sizer_left_panel.Add(lbl_fc, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.text_ctrl_fc, (row, 1), (1, 1), wx.BOTTOM, 5)
        lbl_units_THz = wx.StaticText(self.left_panel, wx.ID_ANY, "(THz):")
        grid_sizer_left_panel.Add(lbl_units_THz, (row, 2), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)

        row += 1
        lbl_laser_bw = wx.StaticText(self.left_panel, wx.ID_ANY, "Laser BW:")
        grid_sizer_left_panel.Add(lbl_laser_bw, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.text_ctrl_laser_bw, (row, 1), (1, 1), wx.BOTTOM, 5)
        lbl_units_laser_bw = wx.StaticText(self.left_panel, wx.ID_ANY, "(x Fc)")
        grid_sizer_left_panel.Add(lbl_units_laser_bw, (row, 2), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)

        row += 1
        lbl_emitted_modes = wx.StaticText(self.left_panel, wx.ID_ANY, "Emitted Modes:")
        lbl_emitted_modes.SetToolTip("The number of emitted modes inside the cavity")
        grid_sizer_left_panel.Add(lbl_emitted_modes, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.text_ctrl_emitted_modes, (row, 1), (1, 1), wx.BOTTOM, 5)

        row += 1
        lbl_reflective_index = wx.StaticText(self.left_panel, wx.ID_ANY, "Reflective Index:")
        grid_sizer_left_panel.Add(lbl_reflective_index, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.text_ctrl_index, (row, 1), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_bandwidth_shape = wx.StaticText(self.left_panel, wx.ID_ANY, "Gain Bandwidth Shape:")
        grid_sizer_left_panel.Add(label_bandwidth_shape, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.combo_bandwidth_shape, (row, 1), (1, 2), wx.BOTTOM, 5)

        row += 1
        grid_sizer_left_panel.Add(self.checkbox_random_phase, (row, 1), (1, 1), wx.LEFT | wx.TOP, 5)

        # SAMPLING PARAMETERS ------------------------------------------------------------------------------------------
        row += 1
        lbl_results = wx.StaticText(self.left_panel, wx.ID_ANY, "SAMPLING")
        lbl_results.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, 0, ""))
        grid_sizer_left_panel.Add(lbl_results, (row, 0), (1, 3), wx.LEFT | wx.RIGHT, 5)

        row += 1
        static_line_3 = wx.StaticLine(self.left_panel, wx.ID_ANY)
        static_line_3.SetMinSize((300, 2))
        grid_sizer_left_panel.Add(static_line_3, (row, 0), (1, 3), wx.BOTTOM | wx.RIGHT | wx.TOP, 5)

        row += 1
        label_window = wx.StaticText(self.left_panel, wx.ID_ANY, "Windowing Function:")
        grid_sizer_left_panel.Add(label_window, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.combo_window, (row, 1), (1, 2), wx.BOTTOM, 5)

        row += 1
        label_Hz = wx.StaticText(self.left_panel, wx.ID_ANY, "Mainlobe Width:")
        grid_sizer_left_panel.Add(label_Hz, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.text_ctrl_mainlobe, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_Hz = wx.StaticText(self.left_panel, wx.ID_ANY, "MLW (Hz)")
        grid_sizer_left_panel.Add(label_Hz, (row, 2), (1, 1), wx.BOTTOM, 5)

        # REPORT -------------------------------------------------------------------------------------------------------
        row += 1

        row += 1
        lbl_results = wx.StaticText(self.left_panel, wx.ID_ANY, "REPORT")
        lbl_results.SetFont(wx.Font(14, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, 0, ""))
        grid_sizer_left_panel.Add(lbl_results, (row, 0), (1, 3), wx.LEFT | wx.RIGHT, 5)

        row += 1
        static_line_3 = wx.StaticLine(self.left_panel, wx.ID_ANY)
        static_line_3.SetMinSize((300, 2))
        grid_sizer_left_panel.Add(static_line_3, (row, 0), (1, 3), wx.BOTTOM | wx.RIGHT | wx.TOP, 5)

        row += 1
        lbl_runtime = wx.StaticText(self.left_panel, wx.ID_ANY, "Total Runtime:")
        grid_sizer_left_panel.Add(lbl_runtime, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_runtime, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_ps = wx.StaticText(self.left_panel, wx.ID_ANY, "(ps)")
        grid_sizer_left_panel.Add(label_ps, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_bandwidth_shape = wx.StaticText(self.left_panel, wx.ID_ANY, "Laser BW:")
        grid_sizer_left_panel.Add(label_bandwidth_shape, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_laser_bw, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_THz = wx.StaticText(self.left_panel, wx.ID_ANY, "(THz)")
        grid_sizer_left_panel.Add(label_THz, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_wavelength = wx.StaticText(self.left_panel, wx.ID_ANY, "Wavelength, λ:")
        grid_sizer_left_panel.Add(label_wavelength, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_wavelength, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_nm = wx.StaticText(self.left_panel, wx.ID_ANY, "(nm)")
        grid_sizer_left_panel.Add(label_nm, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_cavity_modes = wx.StaticText(self.left_panel, wx.ID_ANY, "Cavity Modes, m:")
        grid_sizer_left_panel.Add(label_cavity_modes, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_cavity_modes, (row, 1), (1, 2), wx.BOTTOM, 5)

        row += 1
        label_cavity_length = wx.StaticText(self.left_panel, wx.ID_ANY, "Cavity Length, L:")
        grid_sizer_left_panel.Add(label_cavity_length, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_cavity_length, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_nm = wx.StaticText(self.left_panel, wx.ID_ANY, "(mm)")
        grid_sizer_left_panel.Add(label_nm, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_df = wx.StaticText(self.left_panel, wx.ID_ANY, "Frequency Separation, df:")
        grid_sizer_left_panel.Add(label_df, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_df, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_GHz = wx.StaticText(self.left_panel, wx.ID_ANY, "(GHz)")
        grid_sizer_left_panel.Add(label_GHz, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        label_longitudinal_modes = wx.StaticText(self.left_panel, wx.ID_ANY, "Longitudinal Modes:")
        grid_sizer_left_panel.Add(label_longitudinal_modes, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_longitudinal_modes, (row, 1), (1, 2), wx.BOTTOM, 5)

        row += 1
        lbl_fwhm = wx.StaticText(self.left_panel, wx.ID_ANY, "FWHM:")
        grid_sizer_left_panel.Add(lbl_fwhm, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_fwhm, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_blank = wx.StaticText(self.left_panel, wx.ID_ANY, "")
        grid_sizer_left_panel.Add(label_blank, (row, 2), (1, 1), wx.BOTTOM, 5)

        row += 1
        lbl_fwhm_width = wx.StaticText(self.left_panel, wx.ID_ANY, "FWHM Width:")
        grid_sizer_left_panel.Add(lbl_fwhm_width, (row, 0), (1, 1),
                                  wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.LEFT | wx.RIGHT, 5)
        grid_sizer_left_panel.Add(self.report_fwhm_width, (row, 1), (1, 1), wx.BOTTOM, 5)
        label_GHz = wx.StaticText(self.left_panel, wx.ID_ANY, "(GHz)")
        grid_sizer_left_panel.Add(label_GHz, (row, 2), (1, 1), wx.BOTTOM, 5)

        self.left_panel.SetSizer(grid_sizer_left_panel)

        # PLOT PANEL ===================================================================================================
        grid_sizer_plot.Add(self.canvas, (0, 0), (1, 1), wx.ALL | wx.EXPAND)
        grid_sizer_plot.Add(self.toolbar, (1, 0), (1, 1), wx.ALL | wx.EXPAND)
        grid_sizer_plot.AddGrowableRow(0)
        grid_sizer_plot.AddGrowableCol(0)
        self.plot_panel.SetSizer(grid_sizer_plot)

        # add to main panel --------------------------------------------------------------------------------------------
        grid_sizer_1.Add(self.left_panel, 0, wx.EXPAND | wx.RIGHT, 5)
        grid_sizer_1.Add(self.plot_panel, 1, wx.EXPAND, 5)
        grid_sizer_1.AddGrowableRow(0)
        grid_sizer_1.AddGrowableCol(1)

        sizer_2.Add(grid_sizer_1, 0, wx.EXPAND, 0)

        self.SetSizer(sizer_2)
        self.Layout()

    def popup_dialog(self, message):
        print(message)
        dial = wx.MessageDialog(None, str(message), 'Error', wx.OK | wx.ICON_ERROR)
        dial.ShowModal()

    def get_values(self):
        fc = to_float(self.text_ctrl_fc.GetValue())
        laser_bw = to_float(self.text_ctrl_laser_bw.GetValue())
        emitted_modes = to_integer(self.text_ctrl_emitted_modes.GetValue())
        refraction_index = to_float(self.text_ctrl_index.GetValue())
        bandwidth_shape = str(self.combo_bandwidth_shape.GetValue()).lower()
        window = str(self.combo_window.GetValue()).lower()
        MLW = to_float(self.text_ctrl_mainlobe.GetValue())
        random_phase = bool(self.checkbox_random_phase.GetValue())

        print(fc, emitted_modes, refraction_index, bandwidth_shape, window, MLW, random_phase)

        return fc, laser_bw, emitted_modes, refraction_index, bandwidth_shape, window, MLW, random_phase

    def update(self, evt):
        try:
            params = self.get_values()
            print(params)
            data, plot_data = worker.worker(params)
            self.results_update(data)

            fc, laser_bw, emitted_modes, refraction_index, bandwidth_shape, window, MLW, random_phase = params
            wavelength, laser_bw, df_max, cavity_modes, cavity_length, cavity_df, longitudinal_modes, fwhm_val, fwhm_width, runtime = data

            xt1_left = 0  # show the start of the data
            xt1_right = 2 / cavity_df  # want to display two periods of the frequency step

            xt2_left = 0
            xt2_right = 6 / cavity_df  # what to display two periods of the frequency step

            fc = fc*1e12
            fs = fc * 100

            bound = 2 * cavity_df * (emitted_modes - 1) / 2
            xf_left = max(0.0, fc - bound)
            xf_right = min(fc + bound, fs / 2)

            plot_limits = (xt1_left, xt1_right, xt2_left, xt2_right, xf_left, xf_right)

            self.plot(plot_data, plot_limits)

        except ValueError as e:
            self.popup_dialog(e)

    # ------------------------------------------------------------------------------------------------------------------
    def __do_plot_layout(self):
        self.ax1.set_title('SAMPLED TIMED SERIES DATA')
        self.ax1.set_xlabel('TIME (ps)')
        self.ax1.set_ylabel('AMPLITUDE')

        self.ax2.set_title('SUMMATION OF ALL MODES/TONES')
        self.ax2.set_xlabel('TIME (ps)')
        self.ax2.set_ylabel('AMPLITUDE')

        self.ax3.set_title('SPECTRAL DATA')
        self.ax3.set_xlabel('FREQUENCY (THz)')
        self.ax3.set_ylabel('MAGNITUDE (V)')

        self.ax3.grid()
        self.figure.align_ylabels([self.ax1, self.ax2, self.ax3])
        self.figure.tight_layout()

    def plot(self, plot_data, plot_limits):

        xt, yt_list, yt, yt_envelope, xf_rfft, yf_rfft, yf_smooth = plot_data
        xt1_left, xt1_right, xt2_left, xt2_right, xf_left, xf_right = plot_limits

        xt_scale = 1e12
        xf_scale = 1e12

        # TEMPORAL -----------------------------------------------------------------------------------------------------
        xt_delta = xt[1] - xt[0]
        yt_limit = int(xt1_right / xt_delta)
        yt_limit2 = int(xt2_right / xt_delta)

        self.ax1.clear()
        self.ax1.plot(xt[:yt_limit]*xt_scale, (yt_list[:, :yt_limit]).T)  # All signals
        self.ax1.set_title('SAMPLED TIMED SERIES DATA')
        self.ax1.set_xlabel('TIME (ps)')
        self.ax1.set_ylabel('AMPLITUDE')

        self.temporal_sum.set_data(xt[:yt_limit2] * xt_scale, yt[:yt_limit2])  # The summation of all signals
        self.temporal_hilbert.set_data(xt[:yt_limit2] * xt_scale, yt_envelope[:yt_limit2])  # The envelope of the summation

        self.ax1.set_xlim(left=xt1_left * xt_scale, right=xt1_right * xt_scale)

        self.ax2.set_xlim(left=xt2_left * xt_scale, right=xt2_right * xt_scale)

        # SPECTRAL -----------------------------------------------------------------------------------------------------
        self.spectral.set_data(xf_rfft / xf_scale, np.abs(yf_rfft))  # The spectral plot of sum
        self.spectral_envelope.set_data(xf_rfft / xf_scale, yf_smooth)  # The spectral plot of sum

        self.ax3.set_xlim(left=xf_left / xf_scale, right=xf_right / xf_scale)

        # REDRAW PLOT --------------------------------------------------------------------------------------------------
        self.plot_redraw()

    def plot_redraw(self):
        try:
            self.ax1.relim()  # recompute the ax.dataLim
            self.ax2.relim()  # recompute the ax.dataLim
            self.ax3.relim()  # recompute the ax.dataLim
        except MemoryError as e:
            # xt_length = len(self.ax1.get_xdata())
            # yt_length = len(self.ax1.get_ydata())
            # print(f'Are the lengths of xt: {xt_length} and yt: {yt_length} mismatched?')
            raise ValueError(str(e))
        self.ax1.margins(x=0)
        self.ax1.autoscale(axis='y')
        self.ax2.autoscale(axis='y')
        self.ax3.autoscale(axis='y')

        # UPDATE PLOT FEATURES -----------------------------------------------------------------------------------------
        self.figure.tight_layout()
        self.toolbar.update()  # Not sure why this is needed - ADS
        self.canvas.draw()
        self.canvas.flush_events()

    def results_update(self, data):
        wavelength, laser_bw, df_max, cavity_modes, cavity_length, cavity_df, longitudinal_modes, fwhm_val, fwhm_width, runtime = data

        self.report_runtime.SetValue(str(round(runtime * 1e12, 3)))
        self.report_laser_bw.SetValue(str(laser_bw / 1e12))
        self.report_cavity_length.SetValue(str(laser_bw / 1e12))
        self.report_wavelength.SetValue(str(round(wavelength * 1e9, 3)))

        self.report_df.SetValue(str(round(df_max / 1e9, 3)))
        self.report_cavity_modes.SetValue(str(cavity_modes))
        self.report_cavity_length.SetValue(str(round(cavity_length * 1e3, 3)))
        self.report_longitudinal_modes.SetValue(str(longitudinal_modes))

        self.report_fwhm.SetValue(str(round(fwhm_val, 2)))
        self.report_fwhm_width.SetValue(str(round(fwhm_width * 1e12, 3)))

        print('total runtime:', round(runtime * 1e12, 3), 'ps')
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
        print('FWHM value:', round(fwhm_val, 2))
        print('FWHM width:', round(fwhm_width * 1e12, 3), 'ps')


# FOR RUNNING INDEPENDENTLY ============================================================================================
class MyDemoFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((1200, 705))
        self.panel = MyDemoPanel(self)

        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle("Laser Mode Locking")

    def __do_layout(self):
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel, 1, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(sizer)
        self.Layout()


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyDemoFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True


if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
