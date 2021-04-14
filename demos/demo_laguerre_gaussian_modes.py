import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb, factorial

# https://teaching.smp.uq.edu.au/scims/optics/lgbeam.html

p = 3  # Degree of LG mode
l = 0  # Order of LG mode

w0 = 2.0  # Beam waist
k = 2 * np.pi / 532.0e-9  # Wavenumber of light

zR = k * w0 ** 2.0 / 2  # Calculate the Rayleigh range

# Setup the cartesian grid for the plot at plane z
z = 0.0
xx, yy = np.meshgrid(np.linspace(-5, 5, num=250), np.linspace(-5, 5, num=250))

# Calculate the cylindrical coordinates
r = np.sqrt(xx ** 2 + yy ** 2)
phi = np.arctan2(yy, xx)

# Gaussian Distribution for TEM00
U00 = 1.0 / (1 + 1j * z / zR) * np.exp(-r ** 2.0 / w0 ** 2 / (1 + 1j * z / zR))
w = w0 * np.sqrt(1.0 + z ** 2 / zR ** 2)
R = np.sqrt(2.0) * r / w

# Lpl from OT toolbox (Nieminen et al., 2004)
Lpl = comb(p + l, p) * np.ones(np.shape(R))  # x = R(r, z).^2
for m in range(1, p + 1):
    Lpl = Lpl + (-1.0) ** m / factorial(m) * comb(p + l, p - m) * R ** (2.0 * m)

U = U00 * R ** l * Lpl * np.exp(1j * l * phi) * np.exp(-1j * (2 * p + l + 1) * np.arctan(z / zR))

plt.figure()
plt.title('Intensity')
plt.pcolor(abs(U) ** 2, cmap='plasma')  # 'jet' or 'plasma' or 'nipy_spectral'
plt.axis('equal')

# plt.figure()
# plt.title('Phase')
# plt.pcolor(np.angle(U)**2)
# plt.axis('equal')

plt.show()
