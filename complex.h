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

struct complex
{
	float r;
	float i;
};
 
complex __operator__mul__ (complex c1, complex c2)
{
    complex ret;

    float s1 = c1.r * c2.r;
    float s2 = c1.i * c2.i;
    float s3 = (c1.r + c1.i) * (c2.r + c2.i);

    ret.r = s1 - s2; ret.i = s3 - s1 - s2;

    return ret;
}
complex __operator__mul__ (float s, complex c)
{
    complex ret;
    ret.r = s*c.r;
    ret.i = s*c.i;
    return ret;
}
complex __operator__mul__ (complex c, float s)
{
    complex ret;
    ret.r = s*c.r;
    ret.i = s*c.i;
    return ret;
}

complex __operator__div__ (complex c1, complex c2)
{
    complex ret;

    float s1 = c1.r * c2.r;
    float s2 = c1.i * c2.i;
    float s3 = (c1.r + c1.i) * (c2.r - c2.i);
    float denom = c2.r*c2.r + c2.i*c2.i;

    ret.r = (s1 + s2) / denom; ret.i = (s3 - s1 + s2) / denom;

    return ret;
}
complex __operator__div__ (float s, complex c)
{
    complex ret = complex(s, 0.0);
    ret = ret / c; 
    return ret;
}
complex __operator__div__ (complex c, float s)
{
    complex ret;
    ret = c * (1.0/s); 
    return ret;
}

complex __operator__add__ (complex c1, complex c2)
{
    return complex(c1.r + c2.r, c1.i + c2.i);
}
complex __operator__add__ (complex c, float s)
{
    return complex(c.r + s, c.i);
}
complex __operator__add__ (float s, complex c)
{
    return complex(c.r + s, c.i);
}

complex __operator__sub__ (complex c1, complex c2)
{
    return complex(c1.r - c2.r, c1.i - c2.i);
}
complex __operator__sub__ (complex c, float s)
{
    return complex(c.r - s, c.i);
}
complex __operator__sub__ (float s, complex c)
{
    return complex(s - c.r, -c.i);
}

complex __operator__neg__ (complex c)
{
    return complex(-c.r, -c.i);
}

complex conjugate(complex c)
{
    return complex(c.r, -c.i);
}

float modulus(complex c)
{
    return sqrt(c.r*c.r + c.i*c.i);
}
float modulus_sqrd(complex c)
{
    return c.r*c.r + c.i*c.i;
}

complex pow(complex c, float p)
{
    float rn = pow(modulus(c), p);
    float theta = atan2(c.i, c.r) * p;
    return complex(rn * cos(theta), rn * sin(theta));
}
complex sqrt(complex c)
{
    return pow(c, 0.5);
}
complex exp(complex c)
{
    complex ret;
    float x = exp(c.r);
    ret.r = x * cos(c.i);
    ret.i = x * sin(c.i);
    return ret;
}