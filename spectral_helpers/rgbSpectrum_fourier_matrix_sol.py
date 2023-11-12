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

# title "Rendering biological iridescences with RGB-based renderers"
# url = https://doi.org/10.1145/1122501.1122506
# doi = 10.1145/1122501.1122506

# The solution for equation (5.6)
# Solves for the values x_1, x_2, x_3

from spectral import *
import numpy as np

def B(l):
    x = (2.0 * pi * (l - LAMBDA_MIN)) / (LAMBDA_MAX - LAMBDA_MIN)
    return np.array([1.0,cos(x),sin(x)])

a1 = np.array([0.,0.,0.])
a2 = np.array([0.,0.,0.])
a3 = np.array([0.,0.,0.])

for i in range(LAMBDA_SAMPLES):
    _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
    _B = B(_lambda)
    
    _w = int(float(round(_lambda) - CIE_MIN) / CIE_STEP)
    k = CIE_D65[_w]/CIE_Y_integral

    a1 += k * _B * CIE_X[_w] * LAMBDA_STEP
    a2 += k * _B * CIE_Y[_w] * LAMBDA_STEP
    a3 += k * _B * CIE_Z[_w] * LAMBDA_STEP

mat = np.matrix([a1, a2, a3])
Mf = np.linalg.inv(mat)
Mf = np.transpose(Mf)

print("Mf = "+str(Mf))