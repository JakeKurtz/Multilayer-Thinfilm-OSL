# ------------------------------------------------------------------------- #
#
#    Copyright (C) 2023 Jake Kurtz
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------- #

from math import *
import matplotlib.pyplot as plt

from spectral import *

wavelengths = [0] * LAMBDA_SAMPLES

X_coords = [0] * LAMBDA_SAMPLES
Y_coords = [0] * LAMBDA_SAMPLES
Z_coords = [0] * LAMBDA_SAMPLES

for i in range(LAMBDA_SAMPLES):

    _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
    xyz = CMF_to_XYZ(_lambda)

    wavelengths[i] = _lambda
    X_coords[i] = xyz[0]
    Y_coords[i] = xyz[1]
    Z_coords[i] = xyz[2]

plt.title("CIE Color Matching Functions")

plt.plot(wavelengths, X_coords, 'red', label=r'$\bar{x}(\lambda)$')
plt.plot(wavelengths, Y_coords, 'lime', label=r'$\bar{y}(\lambda)$')
plt.plot(wavelengths, Z_coords, 'blue', label=r'$\bar{z}(\lambda)$')
plt.xlabel(r'$\lambda/nm$')

plt.legend(loc='upper right')

plt.show()