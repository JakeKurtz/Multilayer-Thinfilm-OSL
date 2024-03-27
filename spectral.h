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
#include "util.h"
#pragma once

#define LAMBDA_SAMPLES 12
#define LAMBDA_MIN 380
#define LAMBDA_MAX 780
#define LAMBDA_STEP float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES)

#define CIE_Y_INTEGRAL 106.85691664936203
#define CIE_Y_D65_INTEGRAL 98.89001078341919
#define SCALE float(LAMBDA_MAX - LAMBDA_MIN) / (float(LAMBDA_SAMPLES) * CIE_Y_D65_INTEGRAL)

float pow2(float x)
{
    return x*x;
}
float gauss(float x, float a, float b, float c)
{
    return a*exp(-pow2((x-b)/c));
}
float sigmoidal(float x,float a,float b,float c,float d) 
{
    return d + (a-d) / (1.0 + pow(x/c, b));
}
float sigmoidal(float x,float a,float b,float c)
{
    return a / (1.0 + exp(-b*(x-c)));
}

vector XYZ_to_sRGB(vector xyz)
{
    vector rgb;
    rgb.x = dot(xyz, vector( 3.2404542, -1.5371385, -0.4985314));
    rgb.y = dot(xyz, vector(-0.9692660,  1.8760108,  0.0415560));
    rgb.z = dot(xyz, vector( 0.0556434, -0.2040259,  1.0572252));
    return rgb;
}
vector sRGB_to_XYZ(vector rgb)
{
    vector xyz;
    xyz.x = dot(rgb, vector(0.4124564, 0.3575761, 0.1804375));
    xyz.y = dot(rgb, vector(0.212672,  0.7151522, 0.0721750));
    xyz.z = dot(rgb, vector(0.0193339, 0.1191920, 0.9503041));
    return xyz;
}

vector CMF_to_XYZ(float l)
{
    vector xyz;
    xyz.x = ( l < 500.0 ) ? gauss(l, 0.3762, 449.0573, 27.0226) : gauss(l, 0.8981, 596.7094, 44.6502) + gauss(l, 0.1214, 549.8165, 23.7826);
    xyz.y = gauss(l, 14.0421, 540.7364, 21.3541) + gauss(l, -13.8573, 540.8008, 21.2394) + gauss(l,  0.8341,  559.4890, 60.2144);
    xyz.z = gauss(l, 43.1393, 449.6270, 21.3703) + gauss(l,  0.4729,  469.7765, 38.1569) + gauss(l, -41.5637, 449.6359, 21.1548);
    return xyz;
}

vector SPEC_to_XYZ(float spec[LAMBDA_SAMPLES], float l[LAMBDA_SAMPLES])
{
    
    vector XYZ = vector(0.0, 0.0, 0.0);;
    for (int i = 0; i < LAMBDA_SAMPLES; i++) {
        vector cmf = CMF_to_XYZ(l[i]);
        XYZ += cmf * spec[i];
    }
    return XYZ * SCALE;
}
vector SPEC_to_sRGB(float spec[LAMBDA_SAMPLES], float l[LAMBDA_SAMPLES])
{
    return XYZ_to_sRGB(SPEC_to_XYZ(spec, l));
}
vector SPEC_to_sRGB(float l)
{
    return XYZ_to_sRGB(CMF_to_XYZ(l));
}

/* Fast RGB-to-Spectrum Conversion by Scott Burns */
/* http://scottburns.us/fast-rgb-to-spectrum-conversion-for-reflectances/ */
float sRGB_to_SPEC(vector rgb, float l)
{
    vector opt;
    opt.x = (l > 560.0) ? sigmoidal(l, 0.0144, 135.0296, 590.4639, 0.9761) : sigmoidal(l, 0.0289, 38.6766, 454.7233, 0.0062);
    opt.y = (l < 545.0) ? sigmoidal(l, 0.0123, 63.3450, 488.7451, 0.9711) : sigmoidal(l, 0.0127, -116.3503, 590.4987, 0.9699);
    opt.z = (l < 500.0) ? sigmoidal(l, 0.9622, -0.1542, 489.6339) : sigmoidal(l, 0.0122, -37.0407, 405.3482, 344.2237);
    return dot(rgb, opt);
}

vector Linear_to_sRGB(vector x)
{
  return mix(x*12.92, 1.055 * pow(x, vector(0.41667)) - 0.055, step(0.0031308, x));
}
vector sRGB_to_Linear(vector x)
{
    return mix(x / 12.92, pow((x+0.055)/1.055,vector(2.4)), step(0.04045, x));
}

void lambda_hero(point p, output float lambda_samples[LAMBDA_SAMPLES])
{
    float lambda_r = float(LAMBDA_MAX - LAMBDA_MIN);
    float lambda_h = rand_range(LAMBDA_MIN, LAMBDA_MAX, p);
    for(int i = 0; i < LAMBDA_SAMPLES; i++) {
        float x = (lambda_h - float(LAMBDA_MIN) + (float(i) / float(LAMBDA_SAMPLES)) * lambda_r);
        lambda_samples[i] = mod(x, lambda_r) + float(LAMBDA_MIN);
    }
}
void lambda_uniform(point p, output float lambda_samples[LAMBDA_SAMPLES])
{
    for (int i = 0; i < LAMBDA_SAMPLES; i++) {
        lambda_samples[i] = rand_range(LAMBDA_MIN, LAMBDA_MAX, p+vector(i));
    }
}
void lambda_fixed(point p, output float lambda_samples[LAMBDA_SAMPLES])
{
    for (int i = 0; i < LAMBDA_SAMPLES; i++) {
        lambda_samples[i] = (LAMBDA_MIN + (i * LAMBDA_STEP));
    }
}

void gen_lambda_samples(point p, output float lambda_samples[LAMBDA_SAMPLES])
{
    lambda_hero(p, lambda_samples);
}