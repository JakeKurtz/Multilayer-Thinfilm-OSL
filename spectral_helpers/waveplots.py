# ------------------------------------------------------------------------- #
#
#    Copyright (C) 2023 Jake Kurtz
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------- #

from math import *
import matplotlib.pyplot as plt
import numpy as np

from spectral import *

def plot_cie_color_matching_funcs():
    wavelengths = [0] * LAMBDA_SAMPLES

    X_coords = [0] * LAMBDA_SAMPLES
    Y_coords = [0] * LAMBDA_SAMPLES
    Z_coords = [0] * LAMBDA_SAMPLES

    Xp_coords = [0] * LAMBDA_SAMPLES
    Yp_coords = [0] * LAMBDA_SAMPLES
    Zp_coords = [0] * LAMBDA_SAMPLES

    D65_coords = [0] * LAMBDA_SAMPLES

    for i in range(LAMBDA_SAMPLES):

        _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
        xyz = CMF_to_XYZ_D65(_lambda)

        wavelengths[i] = _lambda

        D65_coords[i] = D65(_lambda)/CIE_Y_integral

        X_coords[i] = xyz[0] 
        Y_coords[i] = xyz[1]
        Z_coords[i] = xyz[2]

        xyz = CMF_to_XYZ_Guass(_lambda)

        Xp_coords[i] = xyz[0] 
        Yp_coords[i] = xyz[1]
        Zp_coords[i] = xyz[2]


    plt.title("CIE Color Matching Functions")

    plt.plot(wavelengths, X_coords, 'red', label=r'$\bar{x}(\lambda)$')
    #plt.plot(wavelengths, Y_coords, 'lime', label=r'$\bar{y}(\lambda)$')
    #plt.plot(wavelengths, Z_coords, 'blue', label=r'$\bar{z}(\lambda)$')

    plt.plot(wavelengths, Xp_coords, 'red', label=r'$\bar{x}(\lambda)$')
    #plt.plot(wavelengths, Yp_coords, 'lime', label=r'$\bar{y}(\lambda)$')
    #plt.plot(wavelengths, Zp_coords, 'blue', label=r'$\bar{z}(\lambda)$')
    #plt.plot(wavelengths, D65_coords, 'black', label=r'$D65(\lambda)$')

    plt.xlabel(r'$\lambda/nm$')

    plt.legend(loc='upper right')

    plt.show()
def plot_srgb_color_matching_funcs():
    wavelengths = [0] * LAMBDA_SAMPLES

    R_coords = [0] * LAMBDA_SAMPLES
    G_coords = [0] * LAMBDA_SAMPLES
    B_coords = [0] * LAMBDA_SAMPLES

    for i in range(LAMBDA_SAMPLES):

        _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))

        xyz = np.array(CMF_to_XYZ(_lambda)) * (D65(_lambda)/CIE_Y_integral)
        rgb = XYZ_to_RGB(xyz * SCALE)

        wavelengths[i] = _lambda

        R_coords[i] = rgb[0] 
        G_coords[i] = rgb[1]
        B_coords[i] = rgb[2]

    plt.title("sRGB Color Matching Functions")

    plt.plot(wavelengths, R_coords, 'red', label=r'$\bar{r}(\lambda)$')
    plt.plot(wavelengths, G_coords, 'lime', label=r'$\bar{g}(\lambda)$')
    plt.plot(wavelengths, B_coords, 'blue', label=r'$\bar{b}(\lambda)$')

    plt.xlabel(r'$\lambda/nm$')

    plt.legend(loc='upper right')

    plt.show()
def plot_srgb_optimal_rgb_components():
    plt.title("Optimal RGB Components")

    plt.plot(CIE_L, CIE_OPTIMAL_R, 'red', label=r'$\rho^R(\lambda)$')
    plt.plot(CIE_L, CIE_OPTIMAL_G, 'lime', label=r'$\rho^G(\lambda)$')
    plt.plot(CIE_L, CIE_OPTIMAL_B, 'blue', label=r'$\rho^B(\lambda)$')

    plt.xlabel(r'$\lambda/nm$')

    plt.legend(loc='upper right')

    plt.show()
def plot_srgb_optimal_rgb_components_fit():

    color = (1.0, 1.0, 1.0)

    wavelengths = [0] * LAMBDA_SAMPLES

    R_coords = [0] * LAMBDA_SAMPLES
    G_coords = [0] * LAMBDA_SAMPLES
    B_coords = [0] * LAMBDA_SAMPLES

    for i in range(LAMBDA_SAMPLES):

        _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))

        wavelengths[i] = _lambda

        R_coords[i] = rFit_Optimal_Fast(_lambda)
        G_coords[i] = gFit_Optimal_Fast(_lambda)
        B_coords[i] = bFit_Optimal_Fast(_lambda)
    
    plt.title("Optimal RGB Components")

    plt.plot(CIE_L, CIE_OPTIMAL_R, 'black', label=r'$\rho^R(\lambda)$')
    plt.plot(CIE_L, CIE_OPTIMAL_G, 'black', label=r'$\rho^G(\lambda)$')
    plt.plot(CIE_L, CIE_OPTIMAL_B, 'black', label=r'$\rho^B(\lambda)$')

    plt.plot(wavelengths, R_coords, 'red', label=r'${Fit}^R(\lambda)$')
    plt.plot(wavelengths, G_coords, 'lime', label=r'${Fit}^G(\lambda)$')
    plt.plot(wavelengths, B_coords, 'blue', label=r'${Fit}^R(\lambda)$')

    plt.xlabel(r'$\lambda/nm$')

    plt.legend(loc='upper right')

    plt.show()

def blahblah():

    color_0 = (1.0, 0.5, 1.0)
    color_1 = (1.0, 1.0, 0.5)
    color_2 = (0.5, 1.0, 1.0)

    wavelengths = [0] * LAMBDA_SAMPLES

    color_0_coords = [0] * LAMBDA_SAMPLES
    color_1_coords = [0] * LAMBDA_SAMPLES
    color_2_coords = [0] * LAMBDA_SAMPLES

    for i in range(LAMBDA_SAMPLES):

        _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))

        wavelengths[i] = _lambda

        color_0_coords[i] = RGB_to_SPEC(color_0, _lambda)
        color_1_coords[i] = RGB_to_SPEC(color_1, _lambda)
        color_2_coords[i] = RGB_to_SPEC(color_2, _lambda)

    plt.plot(wavelengths, color_0_coords, 'magenta', label=r'$Color 0$')
    plt.plot(wavelengths, color_1_coords, 'black', label=r'$Color 1$')
    plt.plot(wavelengths, color_2_coords, 'purple', label=r'$Color 2$')

    plt.xlabel(r'$\lambda/nm$')

    plt.legend(loc='upper right')

    plt.show()

plot_cie_color_matching_funcs()