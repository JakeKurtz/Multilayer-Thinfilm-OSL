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

struct complex_xyz
{
	vector r;
	vector i;
};

vector B(float lambda)
{
    float x = 1.0 / (LAMBDA_MAX - LAMBDA_MIN);

    vector vec;

    vec.x = 1.0;
    vec.y = cos(M_2PI * (lambda - LAMBDA_MIN) * x);
    vec.z = sin(M_2PI * (lambda - LAMBDA_MIN) * x);

    return vec;
}

vector Mf(vector xyz)
{
    vector vec;
    vec.x =  xyz.x*0.01771023 - xyz.y*0.01157347 + xyz.z*0.00377552;
    vec.y =  xyz.x*0.01275868 - xyz.y*0.01898111 + xyz.z*0.00627120;
    vec.z = -xyz.x*0.02633400 + xyz.y*0.02201364 + xyz.z*0.00310862;
    return vec;
}

float sample_xyz(float lambda, vector xyz)
{
    return max(dot(Mf(xyz), B(lambda)), 1e-6) * CIE_Y_integral;
}

complex sample_xyz(float lambda, complex_xyz xyz)
{
    complex c;
    c.r = sample_xyz(lambda, xyz.r);
    c.i = sample_xyz(lambda, xyz.i);
    return c;
}
