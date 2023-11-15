
import cmath
from spectral import *

from dataclasses import dataclass
import numpy as np

@dataclass
class IOR:
    n: np.array((0,0,0))
    k: np.array((0,0,0))

def sample_IOR(wave, ior):
    n = ior.n
    k = ior.k

    _n = rFit_Optimal(wave) * n[0] + gFit_Optimal(wave) * n[1] + bFit_Optimal(wave) * n[2]
    _k = rFit_Optimal(wave) * k[0] + gFit_Optimal(wave) * k[1] + bFit_Optimal(wave) * k[2]

    return complex(_n, _k)