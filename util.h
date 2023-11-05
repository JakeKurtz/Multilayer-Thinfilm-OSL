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

#include "stdosl.h"
#pragma once

vector max(vector a, float b)
{
    vector out = vector(0.0);
    out.x = max(a.x, b);
    out.y = max(a.y, b);
    out.z = max(a.z, b);
    return out;
}

vector clamp(vector x, float a, float b)
{
    return vector(clamp(x.x, a,b), clamp(x.y, a,b), clamp(x.z, a,b));
}

vector smoothstep(float edge0, float edge1, vector x)
{
    vector t = vector(clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0));
    return t * t * (3.0 - 2.0 * t);
}

float fract(float x)
{
    return x - floor(x);
}

vector fract(vector x)
{
    return x - floor(x);
}

float hash13(point p3)
{
    point p = p3;

	p = fract(p * .1031);
    p += dot(p, point(p.z, p.y, p.x) + 31.32);
    return fract((p.x + p.y) * p.z);
}

vector hash33(point p3)
{
    point p = p3;

	p = fract(p * vector(.1031, .1030, .0973));
    p += dot(p, vector(p.y, p.x, p.z)+33.33);
    return fract((vector(p.x, p.x, p.y) + vector(p.y, p.x, p.x))*vector(p.z, p.y, p.x));
}