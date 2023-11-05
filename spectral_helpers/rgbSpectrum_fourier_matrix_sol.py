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
# solves for the values x_1, x_2, x_3

from spectral import *
import numpy as np

def B(l):
    x = 1.0 / (LAMBDA_MAX - LAMBDA_MIN)

    vec = [0,0,0]

    vec[0] = 1.0
    vec[1] = cos(2.0 * pi * (l - LAMBDA_MIN) * x)
    vec[2] = sin(2.0 * pi * (l - LAMBDA_MIN) * x)

    return vec

a1 = [0.,0.,0.]
a2 = [0.,0.,0.]
a3 = [0.,0.,0.]

for i in range(LAMBDA_SAMPLES):
    _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
    _B = B(_lambda)
    xyz = CMF_to_XYZ(_lambda)

    a1[0] += _B[0] * xyz[0] * LAMBDA_STEP
    a1[1] += _B[1] * xyz[0] * LAMBDA_STEP
    a1[2] += _B[2] * xyz[0] * LAMBDA_STEP

    a2[0] += _B[0] * xyz[1] * LAMBDA_STEP
    a2[1] += _B[1] * xyz[1] * LAMBDA_STEP
    a2[2] += _B[2] * xyz[1] * LAMBDA_STEP

    a3[0] += _B[0] * xyz[2] * LAMBDA_STEP
    a3[1] += _B[1] * xyz[2] * LAMBDA_STEP
    a3[2] += _B[2] * xyz[2] * LAMBDA_STEP

mat = np.matrix([a1, a2, a3])
Mf = np.linalg.inv(mat)
Mf = np.transpose(Mf)

print(Mf)