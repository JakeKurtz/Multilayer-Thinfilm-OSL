/* ------------------------------------------------------------------------- *
*
*    Copyright (C) 2023 Jake Kurtz
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program. If not, see <https://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------- */

#include "spectral.h"
#include "complex.h"
#pragma once

struct IOR
{
    color n;
    color k;
};

complex sample_IOR(float lambda, IOR ior)
{
    color n = ior.n;
    color k = ior.k;

    float _n = rFit_Optimal(lambda) * n.r + gFit_Optimal(lambda) * n.g + bFit_Optimal(lambda) * n.b;
    float _k = rFit_Optimal(lambda) * k.r + gFit_Optimal(lambda) * k.g + bFit_Optimal(lambda) * k.b;

    return complex(_n, _k);
}
