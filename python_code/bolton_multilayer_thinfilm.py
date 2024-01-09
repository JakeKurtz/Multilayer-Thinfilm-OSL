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

# title " On the Origin of the Colour of Labradorite"
# url = https://doi.org/10.1002/pssb.19660180123
# doi = 10.1002/pssb.19660180123

# Implementation of equation (6)

#from spectral import *
from thinfilm import Film, modulus_sqrd

import numpy as np
import math
import cmath

def remap(v, l1, h1, l2, h2):
    return l2 + (v-l1)*(h2-l2) / (h1-l1)

def fresnel_dielectric_cos(cosi, eta):
    c = abs(cosi)
    g = eta * eta - 1.0 + c * c
    result = 0

    if (g > 0.0):
        g = math.sqrt(g)
        A = (g - c) / (g + c)
        B = (c * (g + c) - 1.0) / (c * (g - c) + 1.0)
        result = 0.5 * A * A * (1.0 + B * B)
    else:
        result = 1.0

    return result

def cgf(film, phase):
    return (phase**2  * film.v * .5)

def phase(film, wavelength):
    n, cos_theta, _ = film.sample(wavelength)

    #blah = 1.0/1.56
    #sin_theta_i = blah*blah * (1.0 - (film.cos_theta_0*film.cos_theta_0))
    #cos_theta = math.sqrt(1.0 - sin_theta_i)

    #return (4.0 * math.pi * 1.56 * film.d * cos_theta) / wavelength
    return (4.0 * math.pi * n * film.d * cos_theta.real) / wavelength

def bolton_multilayer_thinfilm(wavelength, films):

    film_a, film_b, air, base = films

    phase_a = phase(film_a, wavelength)
    phase_b = phase(film_b, wavelength)

    alpha_a = cgf(film_a, phase_a)
    alpha_b = cgf(film_b, phase_b)

    phi = cmath.exp(-2.0*(alpha_a + alpha_b))
    psi = cmath.exp(-1.0*(alpha_a + alpha_b))

    x = (cmath.exp(-alpha_a) + cmath.exp(-alpha_b))/(1.0 + psi)
    y = (cmath.exp(-alpha_a) - cmath.exp(-alpha_b))/(1.0 - psi)

    M = (phase_a + phase_b) * .5
    N = (phase_a - phase_b) * .5

    f = fresnel_dielectric_cos(air.cos_theta_0, 1.56)

    numer = 1 - cmath.cos(M)*cmath.cos(N)*x + cmath.sin(M)*cmath.sin(N)*y
    denom = (1 + phi - 2.0 * cmath.cos(2.0*M) * psi)
    maximum = 2.0*(1.0+psi)/(1.0-psi)

    r = (1.0 - phi) * numer / (denom*maximum)
    r = math.sqrt(r.real*r.real + r.imag*r.imag)
    return min(max(r+f,0),1)
     
    
