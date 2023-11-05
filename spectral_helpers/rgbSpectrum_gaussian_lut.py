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

# This is my way of generating the LUT (p.117) for the gaussian construction. It even has a pretty plot :)

from spectral import *
import matplotlib.pyplot as plt

_LAMBDA_SAMPLES = 1000
_LAMBDA_MIN = 380
_LAMBDA_MAX = 650
_LAMBDA_STEP = float(_LAMBDA_MAX - _LAMBDA_MIN) / float(_LAMBDA_SAMPLES)

ANGLE_SAMPLES = 30
ANGLE_MIN = 0
ANGLE_MAX = 265
ANGLE_STEP = float(ANGLE_MAX - ANGLE_MIN) / float(ANGLE_SAMPLES)

def plot_cie():
    x_coords = [0] * LAMBDA_SAMPLES
    y_coords = [0] * LAMBDA_SAMPLES
    for i in range(LAMBDA_SAMPLES):
        _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
        xyz = CMF_to_XYZ(_lambda)
        x_coords[i] = xyz[0] / (xyz[0]+xyz[1]+xyz[2])
        y_coords[i] = xyz[1] / (xyz[0]+xyz[1]+xyz[2])

    plt.plot(x_coords, y_coords, 'black', linestyle='dashed')

def rotate_point(xy, angle):
    vec = [0,0]

    vec[0] = cos(angle)*xy[0] + sin(angle)*xy[1]
    vec[1] = -sin(angle)*xy[0] + cos(angle)*xy[1]

    return vec

def angle_sign_dist(alpha, beta):
    phi = abs(beta - alpha) % 360
    d = 360 - phi if (phi > 180) else phi
    return d

def normalize_angle(angle):
    return (angle if (angle >= 0) else (360 - ((-angle) % 360))) % 360

def find_wavelength(D65_xy, WV, WV_len, angle):
    _min_error = 1e32
    _lambda_fin = 0
    _point = WV
    for i in range(_LAMBDA_SAMPLES):
        _lambda = (_LAMBDA_MIN + (i * _LAMBDA_STEP))

        Q_xy = XYZ_coords(CMF_to_XYZ(_lambda))
        
        WQ_dir = [0,0]
        WQ_dir[0] = Q_xy[0]-D65_xy[0]
        WQ_dir[1] = Q_xy[1]-D65_xy[1]

        WQ = [0,0]
        WQ[0] = D65_xy[0] + (WQ_dir[0])
        WQ[1] = D65_xy[1] + (WQ_dir[1])

        dot = WQ_dir[0]*WV_dir[0] + WQ_dir[1]*WV_dir[1]
        det = WQ_dir[0]*WV_dir[1] - WQ_dir[1]*WV_dir[0]
        _angle = normalize_angle(atan2(det, dot) * (180.0/pi))

        a = abs(angle-_angle)#angle_sign_dist(angle,_angle)
        if (a < _min_error):
            _min_error = a
            _lambda_fin = _lambda
            _point = WQ
            #return _lambda_fin, _point

    return  _lambda_fin, _point#LAMBDA_MAX, XYZ_coords(CMF_to_XYZ(LAMBDA_MAX))

D65_xy = [0.31272, 0.32903]
V_xy =   [0.17556, 0.00529]

WV_dir = [0,0]
WV_dir[0] = V_xy[0]-D65_xy[0]
WV_dir[1] = V_xy[1]-D65_xy[1]

WV = [0,0]
WV[0] = D65_xy[0] + (WV_dir[0])
WV[1] = D65_xy[1] + (WV_dir[1])

WV_len = sqrt(WV_dir[0]*WV_dir[0] + WV_dir[1]*WV_dir[1])

q_lut = [0] * (ANGLE_SAMPLES)

for i in range(ANGLE_SAMPLES):
    _angle = (ANGLE_MIN + (i * ANGLE_STEP))

    P = rotate_point(WV_dir, radians(_angle))

    _lambda, T = find_wavelength(D65_xy, WV, WV_len, _angle)
    plt.plot([D65_xy[0], T[0]], [D65_xy[1], T[1]], color=SPEC_to_RGB2(_lambda), marker = 'o')
    plt.text(T[0] * (1 + 0.01), T[1] * (1 + 0.01) , "{:10.2f}".format(_lambda)+"nm", fontsize=10)

    q_lut[i] = _lambda

print(q_lut)

plt.title("Gaussian Construction LUT")

plt.plot(D65_xy[0], D65_xy[1], color=(0.99979, 1.0, 1.0), marker = 'o', markersize=10)
plot_cie()
plt.show()