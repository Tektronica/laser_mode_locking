import numpy as np

SPEED_OF_LIGHT = 299792458  # speed of light (3e8 m/s)

# Cavity properties ----------------------------------------------------------------------------------------------------
"""
Increasing the cavity length decreases the frequency separation
"""
cavity_length = 0.30  # (m)
df = SPEED_OF_LIGHT / (2 * cavity_length)  # Hertz
round_trip_loss = 0.05
fc = 1e9
Q_factor = (fc / df) * (2 * np.pi / round_trip_loss)  # quality factor
df_GHz = df / 1e9

print('cavity length:', cavity_length * 100, 'cm')
print('\tQ factor:', round(Q_factor, 3))
print('\tfrequency separation:', round(df_GHz, 3), 'GHz')
print('\troundtrip-time:', round((1 / df) * 1e9, 3), 'ns')
print()

# HeNe laser -----------------------------------------------------------------------------------------------------------
laser_bw = 1.5e9  # laser bandwidth of HeNe laser (FWHM linewidth)
longitudinal_modes = int(laser_bw / df)  # the number of modes supported by the laser bandwidth
print('laser bandwidth (HeNe):', round(laser_bw / 1e9, 3), 'GHz')
print('longitudinal modes supported:', longitudinal_modes)
print()

# Nd:YAG crystal laser -------------------------------------------------------------------------------------------------
laser_bw = 210e9  # laser bandwidth of Ti:sapphire laser (FWHM linewidth)
longitudinal_modes = int(laser_bw / df)  # the number of modes supported by the laser bandwidth
print('laser bandwidth (Nd:YAG):', round(laser_bw / 1e9, 3), 'GHz')
print('longitudinal modes supported:', longitudinal_modes)
print()

# Ti:sapphire laser ----------------------------------------------------------------------------------------------------
laser_bw = 128e12  # laser bandwidth of Ti:sapphire laser (FWHM linewidth)
longitudinal_modes = int(laser_bw / df)  # the number of modes supported by the laser bandwidth
print('laser bandwidth (HeNe):', round(laser_bw / 1e12, 3), 'THz')
print('longitudinal modes supported:', longitudinal_modes)
print()

wavelength = 1000e-9
f = SPEED_OF_LIGHT / wavelength
print(f * 1e-12, 'THz')
round_trip_time = 1e-9
cavity_length = round_trip_time * SPEED_OF_LIGHT / (2 * wavelength)/1e6
print(cavity_length)
df = SPEED_OF_LIGHT / (2 * cavity_length)  # Hertz

round_trip_loss = 0.05
Q = f * round_trip_time * (2 * np.pi / round_trip_loss)
Q_factor = (f / df) * (2 * np.pi / round_trip_loss)  # quality factor
print('df', round(df, 3), 'Hz')
print(Q / 1e6, '(1e6)')
print(Q_factor / 1e6, '(1e6)')
