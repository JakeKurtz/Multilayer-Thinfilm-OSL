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

void basis(vector n, output vector xp, output vector yp)
{
  //MBR frizvald but attempts to deal with z == -1
  float k = 1.0/max(1.0 + n.z,0.00001);
  float a =  n.y*k;

  float b =  n.y*a;
  float c = -n.x*a;
    
  xp = vector(n.z+b, c, -n.x);
  yp = vector(c, 1.0-b, -n.y);
}

void ggx_anisotropic(float au, float av, vector wo, vector N, vector LN, vector T, point P, output vector wi, output vector lamella_wi) 
{
  vector hash = hashnoise(P);

  float e0 = hash.x;
  float e1 = hash.y;

  float bar = 0;

  if (e0 > 0.0 && e0 < .25) bar = 0.0;
  else if (e0 > .25 && e0 < .75) bar = M_PI;
  else if (e0 >= .75 && e0 <= 1.0) bar = M_2PI;

  float phi =  atan((av/au) * tan(M_2PI * e0)) + bar;
  float A_phi = pow(cos(phi)/au, 2.0) + pow(sin(phi)/av, 2.0);

  float theta = atan(sqrt(-log(e1) / A_phi));

  vector h = vector(
    sin(theta) * cos(phi),
    cos(theta),
    sin(theta) * sin(phi)
  );

  vector B = cross(N,T);

  vector sample = normalize(T*h.x + N*h.y + B*h.z);

  vector sample_l = sample;

  if (N != LN) {
    vector LB, LT;
    basis(LN, LB, LT);
    sample_l = normalize(LT*h.x + LN*h.y + LB*h.z);
  }

  wi = reflect(-wo, sample);
  lamella_wi = reflect(-wo, sample_l);
}

vector ggx(float a, vector wo, vector N, point P) {
  float a2 = a*a;

  vector hash = hashnoise(P);

  float e0 = hash.x;
  float e1 = hash.y;

  float theta = acos(sqrt((1.0 - e0) / (e0 * (a2 - 1.0) + 1.0)));
  float phi = M_2PI * e1;

  vector h = vector(
    sin(theta) * cos(phi),
    cos(theta),
    sin(theta) * sin(phi)
  );

  normal T, B;
  basis(N, T, B);

  vector sample = normalize(T*h.x + N*h.y + B*h.z);
  vector wi = normalize(reflect(-wo, sample));

  return wi;
}

vector sample_wi(float alpha, vector wo, vector N, point P) 
{
  return ggx(alpha, wo, N, P);
}

void sample_wi(float alpha, float anisotropy, float rotation, vector wo, vector N, vector LN, vector T, point P, output vector wi, output vector lamella_wi) 
{
  float alpha_x = max(alpha * alpha * (1.0 + anisotropy), 0.001);
  float alpha_y = max(alpha * alpha * (1.0 - anisotropy), 0.001);
  ggx_anisotropic(alpha_x, alpha_y, wo, N, LN, T, P, wi, lamella_wi);
}