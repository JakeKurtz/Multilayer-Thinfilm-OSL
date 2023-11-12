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

struct matrix22
{
   float t00; float t01;
   float t10; float t11;
};
 
matrix22 matrix22_identity() 
{
    return matrix22(
        1, 0,
        0, 1
    );
}
matrix22 matrix22_zero() 
{
    return matrix22(
        0, 0,
        0, 0
    );
}

matrix22 inv(matrix22 m)
{
    float denom = (m.t00*m.t11 - m.t01*m.t10);
    float det = (denom == 0.0) ? 0.0 : 1.0 / denom;
    return matrix22(
        det * m.t11, -det * m.t01,
       -det * m.t10,  det * m.t00
    );
}

matrix22 __operator__mul__ (matrix22 m1, matrix22 m2)
{
	return matrix22(
        m1.t00*m2.t00 + m1.t01*m2.t10, m1.t00*m2.t01 + m1.t01*m2.t11,
		m1.t10*m2.t00 + m1.t11*m2.t10, m1.t10*m2.t01 + m1.t11*m2.t11
    );
}
matrix22 __operator__mul__ (float s, matrix22 m)
{
    return matrix22(
        s*m.t00, s*m.t01,
        s*m.t10, s*m.t11
    );
}
matrix22 __operator__mul__ (matrix22 m, float s)
{
    return matrix22(
        s*m.t00, s*m.t01,
        s*m.t10, s*m.t11
    );
}

matrix22 __operator__div__ (matrix22 m1, matrix22 m2)
{
    return m1 * inv(m2);
}
matrix22 __operator__div__ (float s, matrix22 m)
{
    return matrix22(
        s/m.t00, s/m.t01,
        s/m.t10, s/m.t11
    );
}
matrix22 __operator__div__ (matrix22 m, float s)
{
    return matrix22(
        m.t00/s, m.t01/s,
        m.t10/s, m.t11/s
    );
}

matrix22 __operator__add__ (matrix22 m1, matrix22 m2)
{
    return matrix22(
        m1.t00+m2.t00, m1.t01+m2.t01,
        m1.t10+m2.t10, m1.t11+m2.t11
    );
}
matrix22 __operator__add__ (matrix22 m, float s)
{
    return matrix22(
        m.t00+s, m.t01+s,
        m.t10+s, m.t11+s
    );
}
matrix22 __operator__add__ (float s, matrix22 m)
{
    return matrix22(
        m.t00+s, m.t01+s,
        m.t10+s, m.t11+s
    );
}

matrix22 __operator__sub__ (matrix22 m1, matrix22 m2)
{
    return matrix22(
        m1.t00-m2.t00, m1.t01-m2.t01,
        m1.t10-m2.t10, m1.t11-m2.t11
    );
}
matrix22 __operator__sub__ (matrix22 m, float s)
{
    return matrix22(
        m.t00-s, m.t01-s,
        m.t10-s, m.t11-s
    );
}
matrix22 __operator__sub__ (float s, matrix22 m)
{
    return matrix22(
        s-m.t00, s-m.t01,
        s-m.t10, s-m.t11
    );
}