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

#define RGB_FIT_FAST

#define LAMBDA_SAMPLES 10
#define LAMBDA_MIN 380
#define LAMBDA_MAX 780
#define LAMBDA_STEP float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES)

#define CIE_Y_integral 106.85691710117203
#define CONSTANT_K 98.8900310736909

#define D65     vector(0.31272, 0.32903, 0.0)
#define VIOLET  vector(0.17556, 0.00529, 0.0)

#define SCALE float(LAMBDA_MAX - LAMBDA_MIN) / float(LAMBDA_SAMPLES * CIE_Y_integral)

#define CIE_D65 vector(0.877739378068311, 0.9234787654138937, 1.0055102172340407)

float gauss(float x, float a, float b, float c)
{
    return a*exp(-pow((x-b)/c,2.0));
}
float sigmoidal(float x,float a,float b,float c,float d) {
    return d + (a-d) / (1.0 + pow(x/c, b));
}
float sigmoidal(float x,float a,float b,float c) {
    return a / (1.0 + exp(-b*(x-c)));
}

float rFit_Optimal(float l) 
{
    #ifdef RGB_FIT_FAST
    return sigmoidal(l, 0.9742, 0.2231, 590.2785);
    #else
    if ( l < 410) {
        return 0.028818291;
    } else if ( l > 410 and  l < 556) {
        float w =  l *0.0219;
        return (0.0224 - 0.0005*cos(w) + 0.0173*sin(w) + 0.0011*cos(2*w) + 0.0056*sin(2*w) - 0.0026*cos(3*w) + 0.0014*sin(3*w));
    } else if ( l > 556 and  l < 584) {
        return gauss( l, 0.0488, 587.7905, 4.2075) + gauss( l, 0.0202, 575.6466, 22.5522) + gauss( l, 0.6490, 606.0866, 17.7174);
    } else if ( l > 584 and  l < 610) {
        return gauss( l, 0.3561, 614.0933, 6.7041) + gauss( l, 0.0950, 593.4713, 5.0340) + gauss( l, 0.8790, 602.8931, 14.7597);
    } else if ( l > 610 and  l < 675) {
        return gauss( l, 0.9752, 681.1453, 374.2223) + gauss( l, -0.0606, 591.3084, 14.8096) + gauss( l, 0.0153, 616.1511, 36.7059);
    } else {
        return 0.976087481;
    }
    #endif
}
float gFit_Optimal(float l)
{
    #ifdef RGB_FIT_FAST
    if (l < 545) return sigmoidal(l, 0.0123, 63.3450, 488.7451, 0.9711);
    else return sigmoidal(l, 0.0127, -116.3503, 590.4987, 0.9699);
    #else
    if ( l < 409) {
        return 0.011877002;
    } else if ( l > 409 and  l < 440) {
        return sigmoidal( l, 0.0118, 56.8665, 444.8670, 0.0196);
    } else if ( l > 440 and  l < 482) {
        return gauss( l, 0.8477, 501.7479, 14.3983) + gauss( l, 235.2388, 483.5420, 0.4566) + gauss( l, 0.0821, 489.6808, 20.0566) + gauss( l, 718.6786, 1.6779e+03, 376.8481);
    } else if ( l > 482 and  l < 514) {
        return gauss( l, 0.9274, 517.2619, 28.0023) + gauss( l, 0.1852, 497.6950, 10.8156) + gauss( l, 0.0858, 490.7573, 7.2268);
    } else if ( l > 514 and  l < 575) {
        return gauss( l, 0.9628, 538.2772, 99.8879) + gauss( l, 0.0299, 508.4227, 17.8645) + gauss( l, 0.0981, 580.3113, 24.7650);
    } else if ( l > 575 and  l < 600) {
        float w =  l * 0.1296;
        return 0.4840 + 0.4195*cos(w) - 0.1922*sin(w) -0.0377*cos(2*w) - 0.0403*sin(2*w) - 0.0079*cos(3*w) + 0.0461*sin(3*w);
    } else if ( l > 600 and  l < 660) {
        return gauss( l, 73.8529, 514.8036, 32.8885) + gauss( l, 0.0067, 609.8295, 14.0107) + gauss( l, 1.8742e+09, -4.3816e+03, 994.0187);
    } else {
        return 0.01268405;
    }
    #endif
}
float bFit_Optimal(float l) 
{
    #ifdef RGB_FIT_FAST
        return sigmoidal(l, 0.9622, -0.1542, 489.6339);
    #else
    if ( l < 400) {
        return 0.959304707;
    } else if ( l > 400 and  l < 440) {
        float w =  l*0.0785;
        return 0.9604 - 0.0012*cos(w) - 0.0006*sin(w) + 0.0001*cos(2*w) + 0.0003*sin(2*w);
    } else if ( l > 440 and  l < 488) {
        return gauss( l, 0.9186, 431.7012, 39.4569) + gauss( l, 0.2310, 479.3085, 13.5087) + gauss( l, 0.4242, 467.5404, 21.5185) + gauss( l, 0.0936, 484.5533, 8.6270 );
    } else if ( l > 488 and  l < 530) {
        return sigmoidal( l, 0.0236, -51.6, 450.733, 30.7708);
    } else if ( l > 530 and  l < 620) {
        float w =  l * 0.0147;
        return 0.1941 + 0.1777*cos(w) -0.1677*sin(w) -0.0008*cos(2*w) -0.0658*sin(2*w);
    } else {
        return 0.011228468;
    }
    #endif
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

vector SPEC_to_RGB2(float l) {
    vector cmf = CMF_to_XYZ(l);
    vector RGB = XYZ_to_RGB(cmf);

    //RGB[0] = min(max(RGB[0],0.0),1.0);
    //RGB[1] = min(max(RGB[1],0.0),1.0);
    //RGB[2] = min(max(RGB[2],0.0),1.0);

    return RGB;
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
    return vector(xyz.x / denom, xyz.y / denom, 0.0);
}
vector XYZ_coords(float l) 
{
    vector xyz = CMF_to_XYZ(l);
    float denom = xyz.x + xyz.y + xyz.z;
    return vector(xyz.x / denom, xyz.y / denom, 0.0);
}