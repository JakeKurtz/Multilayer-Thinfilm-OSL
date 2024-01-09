# ------------------------------------------------------------------------- #
#
#    Copyright (C) 2024 Jake Kurtz
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

import math
import numpy as np
from spectral import *

ANGLE_SAMPLES = 30
ANGLE_MIN = 0
ANGLE_MAX = 265
ANGLE_STEP = float(ANGLE_MAX - ANGLE_MIN) / float(ANGLE_SAMPLES)

D65_xy = np.array([0.31272, 0.32903])
V_xy =   np.array([0.17556, 0.00529])

def mix(v0, v1, t):
  return (1.0 - t) * v0 + t * v1

def clamp(x, minVal, maxVal):
     return min(max(x, minVal), maxVal)

def normalize_angle(angle):
    return (angle if (angle >= 0) else (360 - ((-angle) % 360))) % 360

def S_Gaussain(lamb, lamb_0, chi, xyz):

    w = mix(150.0, 5.0, chi)
    spec = exp(-pow(lamb - lamb_0, 2.0) / (2.0*w*w))

    area_g = 1.0 * math.sqrt(2.0 * math.pi * w*w)

    _xyz = CMF_to_XYZ(lamb_0)

    fuck = xyz[0]+xyz[1]+xyz[2]
    you = _xyz[0]+_xyz[1]+_xyz[2]

    h = fuck/you
    h2 = area_g

    print(spec*h2)

    #print(spec*area_g)
    #print(2948.2862894257373*7)
    
    #0.016779689085165528
    #145855.63591166318

    #print(spec)
    #print(CMF_to_XYZ(lamb_0))
    #print(h)

    #h = 58187.98/w
    #print(h)

    #print(lamb_0, 58187.98/lamb_0)

    #_xyz = CMF_to_XYZ(spec)

    #print((xyz[0]*xyz[1]*xyz[2]))

    #h = 1900.0/(xyz[0]+xyz[1]+xyz[2])
    #h = lamb_0/(xyz[0]*xyz[1]*xyz[2])

    #print(h)

    return spec

def sample_xyz(lamb, xyz):
    
    q_lut = [
        380.0, 465.6, 474.77, 479.63, 483.68, 485.57, 487.73, 489.62, 491.51, 493.67, 495.56, 498.53, 501.5, 505.55, 511.76, 520.6700000000001, 533.63, 544.7, 552.53, 559.55, 564.6800000000001, 568.73, 572.51, 576.56, 581.69, 585.74, 591.6800000000001, 598.7, 611.6600000000001, 649.73
    ]

    P = np.array(XYZ_coords(xyz))
    WP = P - D65_xy
    WV = V_xy - D65_xy

    dot = WP[0]*WV[0] + WP[1]*WV[1]
    det = WP[0]*WV[1] - WP[1]*WV[0]
    angle = normalize_angle(atan2(det, dot) * (180.0/math.pi))

    x = (angle / float(ANGLE_STEP))

    if (angle < ANGLE_MAX):
        min_index = clamp(math.floor(x), 0, 29)
        max_index = clamp(math.ceil(x), 0, 29)
        t = x - math.floor(x)
        
        lamb_0 = mix(q_lut[min_index], q_lut[max_index], t)
        Q = XYZ_coords_lamb(lamb_0)

        chi = np.linalg.norm(D65_xy-P) / np.linalg.norm(D65_xy-Q)

        #if (chi > 0.8):
        return S_Gaussain(lamb, lamb_0, chi, xyz)
    
    return 0
    #return S_Fourier(lamb, xyz)

color_t = np.array([1.0,0.0,0.0])
#sample_xyz(650.0, RGB_to_XYZ(color_t))

#color_t = np.array([0.0,1.0,0.0])

sample_xyz(650.0, RGB_to_XYZ(color_t))
#sample_xyz(532.0, RGB_to_XYZ(color_t))
#sample_xyz(450.0, RGB_to_XYZ(color_t))