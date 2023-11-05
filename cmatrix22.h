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
#include "complex.h"
#pragma once

struct cmatrix22
{
   complex t00; 
   complex t01;
   complex t10; 
   complex t11;
};
 
cmatrix22 identity() 
{
    return cmatrix22(
        complex(1.0, 0.0), complex(0.0, 0.0),
        complex(0.0, 0.0), complex(1.0, 0.0)
    );
}
cmatrix22 zero() 
{
    return cmatrix22(
        complex(0.0, 0.0), complex(0.0, 0.0),
        complex(0.0, 0.0), complex(0.0, 0.0)
    );
}
cmatrix22 one() 
{
    return cmatrix22(
        complex(1.0, 0.0), complex(1.0, 0.0),
        complex(1.0, 0.0), complex(1.0, 0.0)
    );
}

complex det(cmatrix22 m)
{
    return m.t00*m.t11 - m.t01*m.t10;
}
cmatrix22 inv(cmatrix22 m)
{
    complex d = det(m);

    complex det = complex(0.0, 0.0);
    if (modulus_sqrd(d) != 0.0) det = 1.0 / d;

    return cmatrix22(
        det * m.t11, -det * m.t01,
       -det * m.t10,  det * m.t00
    );
}

cmatrix22 __operator__mul__ (cmatrix22 m1, cmatrix22 m2)
{
	return cmatrix22(
        m1.t00*m2.t00 + m1.t01*m2.t10, m1.t00*m2.t01 + m1.t01*m2.t11,
		m1.t10*m2.t00 + m1.t11*m2.t10, m1.t10*m2.t01 + m1.t11*m2.t11
    );
}
cmatrix22 __operator__mul__ (complex s, cmatrix22 m)
{
    return cmatrix22(
        s*m.t00, s*m.t01,
        s*m.t10, s*m.t11
    );
}
cmatrix22 __operator__mul__ (cmatrix22 m, complex s)
{
    return cmatrix22(
        s*m.t00, s*m.t01,
        s*m.t10, s*m.t11
    );
}

cmatrix22 __operator__div__ (cmatrix22 m1, cmatrix22 m2)
{
    return m1 * inv(m2);
}
cmatrix22 __operator__div__ (complex s, cmatrix22 m)
{
    return cmatrix22(
        s/m.t00, s/m.t01,
        s/m.t10, s/m.t11
    );
}
cmatrix22 __operator__div__ (cmatrix22 m, complex s)
{
    return cmatrix22(
        m.t00/s, m.t01/s,
        m.t10/s, m.t11/s
    );
}

cmatrix22 __operator__add__ (cmatrix22 m1, cmatrix22 m2)
{
    return cmatrix22(
        m1.t00+m2.t00, m1.t01+m2.t01,
        m1.t10+m2.t10, m1.t11+m2.t11
    );
}
cmatrix22 __operator__add__ (cmatrix22 m, complex s)
{
    return cmatrix22(
        m.t00+s, m.t01+s,
        m.t10+s, m.t11+s
    );
}
cmatrix22 __operator__add__ (complex s, cmatrix22 m)
{
    return cmatrix22(
        m.t00+s, m.t01+s,
        m.t10+s, m.t11+s
    );
}

cmatrix22 __operator__sub__ (cmatrix22 m1, cmatrix22 m2)
{
    return cmatrix22(
        m1.t00-m2.t00, m1.t01-m2.t01,
        m1.t10-m2.t10, m1.t11-m2.t11
    );
}
cmatrix22 __operator__sub__ (cmatrix22 m, complex s)
{
    return cmatrix22(
        m.t00-s, m.t01-s,
        m.t10-s, m.t11-s
    );
}
cmatrix22 __operator__sub__ (complex s, cmatrix22 m)
{
    return cmatrix22(
        s-m.t00, s-m.t01,
        s-m.t10, s-m.t11
    );
}

cmatrix22 pow(cmatrix22 m, float p)
{
    cmatrix22 ret;
    ret.t00 = pow(m.t00, p); ret.t01 = pow(m.t01, p);
    ret.t10 = pow(m.t10, p); ret.t11 = pow(m.t11, p);
    return ret;
}

void eigen_val(cmatrix22 m, output complex l1, output complex l2)
{
    complex T = m.t00 + m.t11; 
    complex D = m.t00*m.t11 - m.t01*m.t10;

    complex x = T;
    complex y = sqrt(pow(T, 2.0) - 4.0*D);

    l1 = (x + y) * 0.5;
    l2 = (x - y) * 0.5;
}
void eigen_vec(cmatrix22 m, complex l1, complex l2, output complex e1[2], output complex e2[2])
{
    float b = modulus_sqrd(m.t01);
    float c = modulus_sqrd(m.t10);

    if (c != 0.0) {
        e1[0] = l1 - m.t11; e2[0] = l2 - m.t11;
        e1[1] = m.t10;      e2[1] = m.t10;
    } else if (b != 0.0) {
        e1[0] = m.t01;      e2[0] = m.t01;
        e1[1] = l1 - m.t00; e2[1] = l2 - m.t00;
    } else if (b == 0.0 && c == 0.0) {
        e1[0] = complex(1.0, 0.0); e2[0] = complex(0.0, 0.0);
        e1[1] = complex(0.0, 0.0); e2[1] = complex(1.0, 0.0);
    }
}
void diag(cmatrix22 m, output cmatrix22 D, output cmatrix22 P)
{
    complex l1, l2;
    eigen_val(m, l1, l2);   

    complex e1[2]; complex e2[2];
    eigen_vec(m, l1, l2, e1, e2);

    P.t00 = e1[0]; P.t01 = e2[0];
    P.t10 = e1[1]; P.t11 = e2[1];

    D.t00 = l1; D.t01 = complex(0,0);
    D.t10 = complex(0,0); D.t11 = l2;
}