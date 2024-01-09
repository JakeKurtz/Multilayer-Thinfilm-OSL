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

# Helper functions for building spectral shaders.

from cie_colour_matching import *
from math import *
import numpy as np

LAMBDA_SAMPLES = 471
LAMBDA_MIN = 360
LAMBDA_MAX = 830

LAMBDA_STEP = float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES)

SCALE = float(LAMBDA_MAX - LAMBDA_MIN) / (float(LAMBDA_SAMPLES) * 98.89001062030664)

def XYZ_to_sRGB(xyz):
    r =  xyz[0]*3.240479 - xyz[1]*1.537150 - xyz[2]*0.498535
    g = -xyz[0]*0.969256 + xyz[1]*1.875991 + xyz[2]*0.041556
    b =  xyz[0]*0.055648 - xyz[1]*0.204043 + xyz[2]*1.057311
    return [r, g, b]
def sRGB_to_XYZ(rgb):
    x =  rgb[0]*0.412453 + rgb[1]*0.35758 + rgb[2]*0.180423
    y =  rgb[0]*0.212671 + rgb[1]*0.71516 + rgb[2]*0.072169
    z =  rgb[0]*0.019334 + rgb[1]*0.11919 + rgb[2]*0.950227
    return [x, y, z]

def sRGB_to_SPEC(c, w):
    #return c[0] * rFit_Optimal(w) + c[1] * gFit_Optimal(w) + c[2] * bFit_Optimal(w)
    _w = int(float(round(w) - CIE_MIN) / CIE_STEP)
    return c[0] * CIE_OPTIMAL_R[_w] + c[1] * CIE_OPTIMAL_G[_w] + c[2] * CIE_OPTIMAL_B[_w]

def SPEC_to_sRGB(spec, lambda_samples):
    XYZ = [0.0, 0.0, 0.0]
    for i in range(LAMBDA_SAMPLES):
        cmf = CMF_to_XYZ(lambda_samples[i])
        XYZ[0] += cmf[0] * spec[i]
        XYZ[1] += cmf[1] * spec[i]
        XYZ[2] += cmf[2] * spec[i]
    
    XYZ[0] = XYZ[0] * SCALE
    XYZ[1] = XYZ[1] * SCALE
    XYZ[2] = XYZ[2] * SCALE 

    RGB = XYZ_to_sRGB(XYZ)
    RGB[0] = min(max(RGB[0],0.0),1.0)
    RGB[1] = min(max(RGB[1],0.0),1.0)
    RGB[2] = min(max(RGB[2],0.0),1.0)
    return RGB
def SPEC_to_sRGB_lamb(l):
    cmf = CMF_to_XYZ(l)
    RGB = XYZ_to_sRGB(cmf)

    RGB[0] = min(max(RGB[0],0.0),1.0)
    RGB[1] = min(max(RGB[1],0.0),1.0)
    RGB[2] = min(max(RGB[2],0.0),1.0)

    return RGB

def XYZ_coords(xyz):
    denom = xyz[0] + xyz[1] + xyz[2]
    return (xyz[0] / denom, xyz[1] / denom)
def XYZ_coords_lamb(l):
    xyz = CMF_to_XYZ(l)
    denom = xyz[0] + xyz[1] + xyz[2]
    return (xyz[0] / denom, xyz[1] / denom)

def gen_lambda_samples():
    lambda_samples = [0] * LAMBDA_SAMPLES
    for i in range(LAMBDA_SAMPLES):
        lambda_samples[i] = (LAMBDA_MIN + (i * LAMBDA_STEP))
    return lambda_samples
    '''
    lambda_samples = [0] * LAMBDA_SAMPLES
    lambda_r = float(LAMBDA_MAX - LAMBDA_MIN)
    lambda_h = np.random.uniform(low=LAMBDA_MIN, high=LAMBDA_MAX)
    for i in range(LAMBDA_SAMPLES):
        x = (lambda_h - float(LAMBDA_MIN) + (float(i) / float(LAMBDA_SAMPLES)) * lambda_r)
        lambda_samples[i] = np.mod(x, lambda_r) + float(LAMBDA_MIN)
    return lambda_samples
    '''