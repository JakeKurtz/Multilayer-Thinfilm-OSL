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
#include "random.h"
#pragma once

#define LAMBDA_SAMPLES 20
#define LAMBDA_MIN 380
#define LAMBDA_MAX 780
#define LAMBDA_STEP float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES)

#define CIE_Y_integral 106.856895
#define SCALE float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES * CIE_Y_integral)

#define CIE_D65 vector(0.877739378068311, 0.9234787654138937, 1.0055102172340407)

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
        lambda_samples[i] = rand_range(LAMBDA_MIN, LAMBDA_MAX, p);
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

float gaussian(float x, float mu, float sigma)
{
    return 1.0 / (sigma * sqrt(2.0 * M_PI)) * exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
}

vector CMF_to_XYZ(float l)
{
	return vector(
    	8233.3080 * gaussian(l, 593.949462640494, 34.00) + 1891.2652 * gaussian(l, 448.8951, 18.7851),
        10522.6505 * gaussian(l, 555.3855, 40.7979),
        11254.7819 * gaussian(l, 452.9834, 21.5712)
    ) / 100.;
}

vector XYZ_to_RGB(vector xyz)
{
    float r =  xyz.x*3.240479 - xyz.y*1.537150 - xyz.z*0.498535;
    float g = -xyz.x*0.969256 + xyz.y*1.875991 + xyz.z*0.041556;
    float b =  xyz.x*0.055648 - xyz.y*0.204043 + xyz.z*1.057311;
    return vector(r, g, b);
}
vector RGB_to_XYZ(vector rgb)
{
    float x =  rgb.x*0.412453 + rgb.y*0.35758 + rgb.z*0.180423;
    float y =  rgb.x*0.212671 + rgb.y*0.71516 + rgb.z*0.072169;
    float z =  rgb.x*0.019334 + rgb.y*0.11919 + rgb.z*0.950227;
    return vector(x, y, z);
}

vector SPEC_to_XYZ(float spec[LAMBDA_SAMPLES], float lamb_samples[LAMBDA_SAMPLES])
{
    vector XYZ = vector(0.0, 0.0, 0.0);
    for (int i = 0; i < LAMBDA_SAMPLES; i++) {
        vector cmf = CMF_to_XYZ(lamb_samples[i]);
        XYZ.x += cmf.x * spec[i];
        XYZ.y += cmf.y * spec[i];
        XYZ.z += cmf.z * spec[i];
    }
    return XYZ * SCALE;
}
vector SPEC_to_RGB(float spec[LAMBDA_SAMPLES], float lamb_samples[LAMBDA_SAMPLES])
{
    vector XYZ = vector(0.0, 0.0, 0.0);
    for (int i = 0; i < LAMBDA_SAMPLES; i++) {
        vector cmf = CMF_to_XYZ(lamb_samples[i]);
        XYZ.x += cmf.x * spec[i];
        XYZ.y += cmf.y * spec[i];
        XYZ.z += cmf.z * spec[i];
    }
    vector RGB = XYZ_to_RGB(XYZ * SCALE);
    return RGB;
}
vector SPEC_to_RGB(float lambda)
{
    vector XYZ = vector(0.0, 0.0, 0.0);
    vector cmf = CMF_to_XYZ(lambda);
    vector RGB = XYZ_to_RGB(cmf);
    return RGB;
}

vector XYZ_coords(vector xyz) 
{
    float denom = xyz.x + xyz.y + xyz.z;
    return vector(clamp(xyz.x / denom, 0.172, 0.8), xyz.y / denom, 0.0);
}
vector XYZ_coords(float l) 
{
    vector xyz = CMF_to_XYZ(l);
    float denom = xyz.x + xyz.y + xyz.z;
    return vector(clamp(xyz.x / denom, 0.172, 0.8), xyz.y / denom, 0.0);
}