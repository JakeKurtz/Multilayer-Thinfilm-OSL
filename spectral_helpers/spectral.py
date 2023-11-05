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

LAMBDA_SAMPLES = 1000
LAMBDA_MIN = 380
LAMBDA_MAX = 780
LAMBDA_STEP = float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES)

CIE_Y_integral = 106.856895
SCALE = float(LAMBDA_MAX - LAMBDA_MIN) / (float(LAMBDA_SAMPLES) * CIE_Y_integral)

def XYZ_to_RGB(xyz):
    r =  xyz[0]*3.240479 - xyz[1]*1.537150 - xyz[2]*0.498535
    g = -xyz[0]*0.969256 + xyz[1]*1.875991 + xyz[2]*0.041556
    b =  xyz[0]*0.055648 - xyz[1]*0.204043 + xyz[2]*1.057311
    return [r, g, b]
def RGB_to_XYZ(rgb):
    x =  rgb[0]*0.412453 + rgb[1]*0.35758 + rgb[2]*0.180423
    y =  rgb[0]*0.212671 + rgb[1]*0.71516 + rgb[2]*0.072169
    z =  rgb[0]*0.019334 + rgb[1]*0.11919 + rgb[2]*0.950227
    return [x, y, z]

def SPEC_to_RGB(spec, lambda_samples):
    XYZ = [0.0, 0.0, 0.0]
    for i in range(LAMBDA_SAMPLES):
        cmf = CMF_to_XYZ(lambda_samples[i])
        XYZ[0] += cmf[0] * spec[i]
        XYZ[1] += cmf[1] * spec[i]
        XYZ[2] += cmf[2] * spec[i]
    
    XYZ[0] = XYZ[0] * SCALE
    XYZ[1] = XYZ[1] * SCALE
    XYZ[2] = XYZ[2] * SCALE 

    RGB = XYZ_to_RGB(XYZ)
    RGB[0] = min(max(RGB[0],0.0),1.0)
    RGB[1] = min(max(RGB[1],0.0),1.0)
    RGB[2] = min(max(RGB[2],0.0),1.0)
    return RGB
def SPEC_to_RGB2(l):
    cmf = CMF_to_XYZ(l)
    RGB = XYZ_to_RGB(cmf)

    RGB[0] = min(max(RGB[0],0.0),1.0)
    RGB[1] = min(max(RGB[1],0.0),1.0)
    RGB[2] = min(max(RGB[2],0.0),1.0)

    return RGB

def XYZ_coords(xyz):
    denom = xyz[0] + xyz[1] + xyz[2]
    return (xyz[0] / denom, xyz[1] / denom)

print(XYZ_to_RGB([0.877739378068311, 0.9234787654138937, 1.0055102172340407]))