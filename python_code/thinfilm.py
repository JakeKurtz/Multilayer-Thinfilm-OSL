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


import cmath
import math
import cmath
import numpy as np

from ior import IOR, sample_IOR
from spectral import *
from box_muller import rand

from dataclasses import dataclass

@dataclass
class Film:
    d: float # Thickness (nm)
    v: float # Variance (nm)
    ior: IOR # Index of Refraction

    cos_theta_0: complex
    ior_0: IOR

    def sample(self, wavelength):
        n = sample_IOR(wavelength, self.ior)
        n_0 = sample_IOR(wavelength, self.ior_0)

        cos_theta = cos_theta_i(n_0, n, self.cos_theta_0)
        d = d_i(self.d, self.v)

        return n, cos_theta, d

def conjugate(c):
    return complex(c.real, -c.imag)

def modulus_sqrd(c):
    #return cmath.sqrt(conjugate(c) * c)
    return c.real*c.real + c.imag*c.imag

# The Reflectance for p-polarized light #
def rp(n_m, n_l, cos_theta_m, cos_theta_l):
    return ( 
        n_l*cos_theta_m - n_m*cos_theta_l ) / ( 
        n_m*cos_theta_l + n_l*cos_theta_m )

# The Reflectance for s-polarized light #
def rs(n_m, n_l, cos_theta_m, cos_theta_l):
    return (
        n_m*cos_theta_m - n_l*cos_theta_l) / (
        n_m*cos_theta_m + n_l*cos_theta_l )

# The Transmittance for p-polarized light #
def tp(n_m, n_l, cos_theta_m, cos_theta_l):
    return 2.0*n_m*cos_theta_m / (
        n_m*cos_theta_l + n_l*cos_theta_m)

# The Transmittance for s-polarized light #
def ts(n_m, n_l, cos_theta_m, cos_theta_l):
    return 2.0*n_m*cos_theta_m / (
        n_m*cos_theta_m + n_l*cos_theta_l)

def fresnel(n_m, n_l, cos_theta_m, cos_theta_l,):
    a = n_m*cos_theta_m
    b = n_m*cos_theta_l
    c = n_l*cos_theta_l
    d = n_l*cos_theta_m

    denom_1 = b + d
    denom_2 = a + c

    numer_1 = 2.0*a

    rp = (d - b) / denom_1
    rs = (a - c) / denom_2

    tp = numer_1 / denom_1
    ts = numer_1 / denom_2

    return rp, rs, tp, ts

def cos_theta_i(n_0, n_i, cos_theta_0):
    blah = n_0/n_i
    sin_theta_i = blah*blah * (1.0 - (cos_theta_0*cos_theta_0))
    return cmath.sqrt(1.0 - sin_theta_i)

def d_i(d, v):
    return rand(d, v)

def D_mat(r_ij, r_ji, t_ij, t_ji):
    return (1.0 / t_ij) * np.matrix(
           [[complex(1.0,  0.0), -r_ji],
           [r_ij, (t_ij*t_ji - r_ij*r_ji)]]
    )

def P_mat(wave, ior, d, cos_theta):
    phi = ior * complex(0.0, cos_theta.real * (2.0 * math.pi/wave) * d)

    return np.matrix(
        [[cmath.exp(-phi), complex(0.0, 0.0)],
        [complex(0.0, 0.0), cmath.exp(phi)]]
    )