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

from ior import *
from thinfilm import *

import matplotlib.pyplot as plt

import numpy as np

def compute_interface(wavelength, film_i, film_j):
    
    ior_i, ct_i, d_i = film_i.sample(wavelength)
    ior_j, ct_j, _ = film_j.sample(wavelength)

    rp_ij, rs_ij, tp_ij, ts_ij = fresnel(
        ior_i, ior_j, ct_i, ct_j
    )
    rp_ji, rs_ji, tp_ji, ts_ji = fresnel(
        ior_j, ior_i, ct_j, ct_i
    )

    Ds_ij = D_mat(rs_ij, rs_ji, ts_ij, ts_ji)
    Dp_ij = D_mat(rp_ij, rp_ji, tp_ij, tp_ji)
    P_i = P_mat(wavelength, ior_i, d_i, ct_i)

    return Ds_ij, Dp_ij, P_i

def tmm_multilayer_thinfilm(layers, wavelength, films):

    film_a, film_b, air, base = films

    Ts = np.matrix([
        [complex(1,0), complex(0,0)], 
        [complex(0,0), complex(1,0)]
    ])
    Tp = np.matrix([
        [complex(1,0), complex(0,0)], 
        [complex(0,0), complex(1,0)]
    ])

    # ----------------------------------- Air ---------------------------------- #

    Ds_ij, Dp_ij, P_i = compute_interface(wavelength, air, film_a)
    Ts = Ts*P_i*Ds_ij
    Tp = Tp*P_i*Dp_ij

    # --------------------------------- Layers --------------------------------- #

    for i in range(1,layers-1):
        if (i % 2 == 0):
            Ds_12, Dp_12, P_1 = compute_interface(wavelength, film_a, film_b)
            Ts = Ts*P_1*Ds_12
            Tp = Tp*P_1*Dp_12
        else: 
            Ds_21, Dp_21, P_2 = compute_interface(wavelength, film_b, film_a)
            Ts = Ts*P_2*Ds_21
            Tp = Tp*P_2*Dp_21

    # ---------------------------------- Base ---------------------------------- #
 
    Ds_ij, Dp_ij, P_i = compute_interface(wavelength, film_a, base) if (layers-1 % 2 == 0) else compute_interface(wavelength, film_b, base)

    Ts = Ts*P_i*Ds_ij
    Tp = Tp*P_i*Dp_ij

    # ------------------------------- Send it B) ------------------------------- #

    r_s = Ts[1,0] / Ts[0,0]
    r_p = Tp[1,0] / Tp[0,0]

    t_s = 1.0 / Ts[0,0]
    t_p = 1.0 / Tp[0,0]

    R_s = modulus_sqrd(r_s)
    R_p = modulus_sqrd(r_p)

    base_n, base_ct, _ = base.sample(wavelength)
    air_n, _, _ = air.sample(wavelength)

    x = base_n * base_ct
    y = conjugate(base_n) * base_ct
    z = air_n.real * air.cos_theta_0.real

    T_s = (modulus_sqrd(t_s) * x.real) / z
    T_p = (modulus_sqrd(t_p) * y.real) / z

    R = (R_s + R_p) * 0.5
    T = (T_s + T_p) * 0.5

    return R

