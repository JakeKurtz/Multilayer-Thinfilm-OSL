/* ------------------------------------------------------------------------- *
*
*    Copyright (C) 2024 Jake Kurtz
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

#define INFINITE 1e31
#define EPSILON 1e-6

/* ---------------------------------- Math ---------------------------------- */

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

/* --------------------------------- Random --------------------------------- */

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

void rng_seed(output int rng, int seed)
{
  int chash = seed;
  if (chash == 0) chash = 1;
  rng = chash * 891694213;
}

float rng_uniform(output int rng)
{
  float res = rng / float(2137483647) * 0.5 + 0.5;
  rng *= 891694213;
  return res;
}

float rand_range(float start, float end, vector p)
{
  vector hash = hashnoise(p);
  return (end-start)*hash.x + start;
}

/* Box Muller Transform */
float rand(float mu, float var, vector p)
{
  float u1 = hashnoise(p);
  float u2 = hashnoise(p+vector(u1));
  float r = sqrt(-2.0 * log(u2));
  float t = M_2PI * u2;

  float x = r * cos(t);

  return sqrt(var) * x + mu;
}

/* ----------------------------------- 3D ----------------------------------- */

/* SPDX-FileCopyrightText: 2009-2010 Sony Pictures Imageworks Inc., et al. All Rights Reserved.
 * SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Adapted code from Open Shading Language. */
float fresnel(float cosi, float eta)
{
  /* compute fresnel reflectance without explicitly computing
   * the refracted direction */
  float c = fabs(cosi);
  float g = eta * eta - 1 + c * c;
  float result;

  if (g > 0.0) {
    g = sqrt(g);
    float A = (g - c) / (g + c);
    float B = (c * (g + c) - 1) / (c * (g - c) + 1);
    result = 0.5 * A * A * (1 + B * B);
  }
  else {
    result = 1.0; /* TIR (no refracted component) */
  }

  return result;
}

/* Rotate a vector perpendicular to another */
/* https://math.stackexchange.com/questions/3130813/rotating-a-vector-perpendicular-to-another */
vector rotate_perp(point Q, float angle, vector axis)
{
  vector z = cross(axis, Q);
  return cos(angle) * Q + sin(angle) * z;
}

/* https://www.shadertoy.com/view/ltdXRN */
void calc_tangent_space( vector N, float rotation, output vector TangentU, output vector TangentV )
{
	vector Up = vector( 0., 0., 1. );
  TangentV = normalize(rotate_perp(cross(Up, N), rotation, N));
  TangentU = cross( TangentV, N );
}

/* Sampling Anisotropic Microfacet BRDF */
/* https://agraphicsguynotes.com/posts/sample_anisotropic_microfacet_brdf/ */
void ggx_anisotropic(float au, float av, float rotation, vector wo, vector N, vector LN, point P, output vector wi, output vector lamella_wi) 
{
  float e0 = hashnoise(P);
  float e1 = hashnoise(P+wo);

  //vector hash = hashnoise(P);
  //float e0 = hash.x;
  //float e1 = hash.y;

  float bar = 0.0;

  if (e0 > .25 && e0 < .75) bar = M_PI;
  else if (e0 >= .75 && e0 <= 1.0) bar = M_2PI;

  float phi = atan((av/au) * tan(M_2PI * e0)) + bar;

  float cos_phi = cos(phi); float sin_phi = sin(phi);

  float A_phi = pow(cos_phi/au, 2.0) + pow(sin_phi/av, 2.0);

  //float theta = atan(sqrt(e1 / ((1.0 - e1)*A_phi)));
  float theta = acos(sqrt((A_phi*(1.0-e1)) / (e1*(1.0-A_phi)+A_phi)));

  float sin_theta = sin(theta);

  vector h = vector(
    sin_theta * cos_phi,
    cos(theta),
    sin_theta * sin_phi
  );

  normal B, T;
  calc_tangent_space(N, rotation, B, T);

  vector sample = normalize(T*h.x + N*h.y + B*h.z);
  vector sample_l = sample;

  if (N != LN) {
    normal LB, LT;
    calc_tangent_space(LN, rotation, LB, LT);
    sample_l = normalize(LT*h.x + LN*h.y + LB*h.z);
  }

  wi = reflect(-wo, sample);
  lamella_wi = reflect(-wo, sample_l);
}

vector ggx(float a, vector wo, vector N, point P) {
  float a2 = a*a;

  float e0 = hashnoise(P);
  float e1 = hashnoise(P+wo);

  float theta = acos(sqrt((1.0 - e0) / (e0 * (a2 - 1.0) + 1.0)));
  float phi = M_2PI * e1;

  vector h = vector(
    sin(theta) * cos(phi),
    cos(theta),
    sin(theta) * sin(phi)
  );

  normal B, T;
  calc_tangent_space(N, 0.0, B, T);

  vector sample = normalize(T*h.x + N*h.y + B*h.z);
  vector wi = normalize(reflect(-wo, sample));

  return wi;
}

vector sample_wi(float alpha, vector wo, vector N, point P) 
{
  return ggx(alpha, wo, N, P);
}

void sample_wi(float alpha, float anisotropy, float rotation, vector wo, vector N, vector LN, point P, output vector wi, output vector lamella_wi) 
{
  float alpha_x = max(alpha * alpha * (1.0 + anisotropy), 0.001);
  float alpha_y = max(alpha * alpha * (1.0 - anisotropy), 0.001);
  ggx_anisotropic(alpha_x, alpha_y, rotation, wo, N, LN, P, wi, lamella_wi);
}