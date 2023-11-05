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

#include "complex.h"
#include "cmatrix22.h"

#define INFINITE 1e31

/* The Reflectance for p-polarized light */
complex rp(complex n_m, complex n_l, complex cos_theta_m, complex cos_theta_l)
{
    return ( 
        n_l*cos_theta_m - n_m*cos_theta_l ) / ( 
        n_m*cos_theta_l + n_l*cos_theta_m );
}

/* The Reflectance for s-polarized light */
complex rs(complex n_m, complex n_l, complex cos_theta_m, complex cos_theta_l)
{
    return (
        n_m*cos_theta_m - n_l*cos_theta_l) / (
        n_m*cos_theta_m + n_l*cos_theta_l );
}

/* The Transmittance for p-polarized light */
complex tp(complex n_m, complex n_l, complex cos_theta_m, complex cos_theta_l)
{
    return 2.0*n_m*cos_theta_m / (
        n_m*cos_theta_l + n_l*cos_theta_m);
}

/* The Transmittance for s-polarized light */
complex ts(complex n_m, complex n_l, complex cos_theta_m, complex cos_theta_l)
{
    return 2.0*n_m*cos_theta_m / (
        n_m*cos_theta_m + n_l*cos_theta_l);
}

void compute_polarization(
    complex n_m, complex n_l, complex cos_theta_m, complex cos_theta_l,
    output complex rp, 
    output complex rs, 
    output complex ts, 
    output complex tp) 
{
    complex a = n_m*cos_theta_m;
    complex b = n_m*cos_theta_l;
    complex c = n_l*cos_theta_l;
    complex d = n_l*cos_theta_m;

    complex denom_1 = b + d;
    complex denom_2 = a + c;

    complex numer_1 = 2.0*a;

    rp = (d - b) / denom_1;
    rs = (a - c) / denom_2;

    tp = numer_1 / denom_1;
    ts = numer_1 / denom_2;
}

complex cos_theta_i(complex n, complex n_0, float sin_theta_0)
{
    complex sin_theta = (conjugate(n) / pow(modulus(n), 2.0)) * n_0 * sin_theta_0;
    return sqrt(1.0 - pow(sin_theta, 2.0));
}

cmatrix22 D_mat(complex r_ij, complex r_ji, complex t_ij, complex t_ji)
{
    return (1.0 / t_ij) * cmatrix22(
           complex(1.0,  0.0), -r_ji,
           r_ij, (t_ij*t_ji - r_ij*r_ji)
    );
}

cmatrix22 P_mat(float lambda, complex ior, float d, complex cos_theta)
{   
    complex phi;
    if (d >= INFINITE) {
        phi = ior.r * cos_theta.r * complex(0.0, (M_2PI/lambda) * d);
    } else {
        phi = ior * cos_theta * complex(0.0, (M_2PI/lambda) * d);
    }

    return cmatrix22(
        exp(-phi), complex(0.0, 0.0),
        complex(0.0, 0.0), exp(phi)
    );
}