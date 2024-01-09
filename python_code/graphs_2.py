from thinfilm import Film, cos_theta_i
from ior import IOR
import math
import numpy as np
import matplotlib.pyplot as plt

from spectral import *

from tmm_multilayer_thinfilm import tmm_multilayer_thinfilm
from bolton_multilayer_thinfilm import bolton_multilayer_thinfilm

INFINITE = 1e31

def f_sd(wavelength, sd, films_tmm, films_bol):

    R_1 = [0] * len(wavelength)
    R_2 = [0] * len(wavelength)

    for i in range(len(wavelength)):
        w = wavelength[i]

        #films_bol[0].v = sd[i]**2 / films_bol[0].d**2

        R_1[i] = bolton_multilayer_thinfilm(w, films_bol)
        
        for j in range(1):
            R_2[i] += tmm_multilayer_thinfilm(100, w, films_tmm)/1.0

    return R_1, R_2

I = 45.0 # angle of incidence in degrees

ior_0 = IOR(1.0, 0.0)
cos_theta_0 = math.cos(math.radians(I))

a_m = 90.6 # lamellae mean thickness
a_sd = 1.0 # lamellae standard deviation of thickness
a_v = a_sd**2
a_v_bol = 8.0**2.0 / a_m**2

b_m = 60.0 # lamellae mean thickness
b_sd = 1.0 # lamellae standard deviation of thickness
b_v = b_sd**2
b_v_bol = 8.0**2.0 / b_m**2

film_a = Film(a_m, a_v, IOR(1.6, 0.0), cos_theta_0, ior_0)
film_b = Film(b_m, b_v, IOR(1.5, 0.0), cos_theta_0, ior_0)

film_a_bol = Film(a_m, a_v_bol, IOR(1.6, 0.0), cos_theta_0, ior_0)
film_b_bol = Film(b_m, b_v_bol, IOR(1.5, 0.0), cos_theta_0, ior_0)

air = Film(INFINITE, 1.0, ior_0, cos_theta_0, ior_0)
base = Film(INFINITE, 1.0, IOR(1.5, 0.0), cos_theta_0, ior_0)

films_bol = [film_a_bol, film_b_bol, air, base]
films_tmm = [film_a, film_b, air, base]

wavelength = CIE_L

min_angle = 0.0
max_angle = 90.0
angle_step = float(max_angle - min_angle) / float(CIE_SAMPLES-1)
angle = [x * angle_step for x in range(CIE_SAMPLES)]

min_d = 50.0
max_d = 150.0
d_step = float(max_d - min_d) / float(CIE_SAMPLES-1)
thickness = [(x * d_step + min_d) for x in range(CIE_SAMPLES)]

min_sd = 0.00
max_sd = 500.0
sd_step = float(max_sd - min_sd) / float(CIE_SAMPLES-1)
sd = [(x * sd_step + min_sd) for x in range(CIE_SAMPLES)]


X, Y = np.meshgrid(wavelength, [30])
foo, bar = np.array(f_sd(np.ravel(X), np.ravel(Y), films_tmm, films_bol))

plt.plot(wavelength, foo, 'black', label=r'${Bolton}(\lambda)$')
plt.plot(wavelength, bar, 'red', label=r'${TMM}(\lambda)$')
plt.xlabel(r'$\lambda/nm$')
plt.legend(loc='upper right')
plt.show()

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(wavelength, sd)
zs = np.array(f_sd(np.ravel(X), np.ravel(Y), films_tmm, films_bol))

Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel(r'$\lambda/nm$')
ax.set_ylabel(r'$\theta$')
ax.set_zlabel('Z Label')

plt.xlim([CIE_MIN, CIE_MAX])
#ax.set_zlim([0, 1])
plt.gca().invert_yaxis()

plt.show()
'''
