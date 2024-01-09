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

from dataclasses import dataclass
import numpy as np

from spectral import *

@dataclass
class IOR:
    n: np.array((0,0,0))
    k: np.array((0,0,0))

    def __init__(self, _n: float, _k: float):
        self.n = np.array((_n,_n,_n))
        self.k = np.array((_k,_k,_k))

def sample_IOR(wave, ior):
    n = ior.n
    k = ior.k

    _n = rFit_Optimal(wave) * n[0] + gFit_Optimal(wave) * n[1] + bFit_Optimal(wave) * n[2]
    _k = rFit_Optimal(wave) * k[0] + gFit_Optimal(wave) * k[1] + bFit_Optimal(wave) * k[2]
    return complex(_n, _k)