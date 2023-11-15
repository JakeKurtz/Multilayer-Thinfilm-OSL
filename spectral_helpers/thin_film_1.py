from ior import *
from thinfilm import *
import cmath

import matplotlib.pyplot as plt

import numpy as np
from math import *

LAYERS = 3
INFINITE = 1e31

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def normalize(v):
    return v / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def reflect(I, N):
    return I - 2.0 * np.dot(N, I) * N

def basis(n):
    #MBR frizvald but attempts to deal with z == -1
    k = 1.0/max(1.0 + n[2],0.00001)
    a =  n[1]*k

    b =  n[1]*a
    c = -n[0]*a

    xp = np.array((n[2]+b, c, -n[0]))
    yp = np.array((c, 1.0-b, -n[1]))

    return xp, yp

def sample_wi(r, wo, n):
    a2 = r*r*r*r

    e0 = np.random.rand(1,1)
    e1 = np.random.rand(1,1)

    theta = acos(sqrt((1.0 - e0) / (e0 * (a2 - 1.0) + 1.0)))
    phi = 2.0 * pi * e1

    h = np.array((
        sin(theta) * cos(phi),
        cos(theta),
        sin(theta) * sin(phi)
    ))

    t, b = basis(n)

    sample = normalize(t * h[0] + n * h[1] + b * h[2])
    wi = normalize(reflect(-wo, sample))
    
    return wi

def thin_film(wave, n, d, cos_theta):
    Ts = np.matrix([[complex(1,0), complex(0,0)], [complex(0,0), complex(1,0)]])
    Tp = np.matrix([[complex(1,0), complex(0,0)], [complex(0,0), complex(1,0)]])

    n_i = complex(0,0)
    ct_i = complex(0,0)
    n_j = complex(0,0)
    ct_j = complex(0,0)

    rs_ij= complex(0,0) 
    ts_ij= complex(0,0)
    rs_ji= complex(0,0) 
    ts_ji= complex(0,0)
    rp_ij= complex(0,0) 
    tp_ij= complex(0,0)
    rp_ji= complex(0,0) 
    tp_ji= complex(0,0)

    for i in range(LAYERS-1):
        n_i = n[i];   ct_i = cos_theta[i]
        n_j = n[i+1]; ct_j = cos_theta[i+1]

        rp_ij, rs_ij, tp_ij, ts_ij = compute_polarization(n_i, n_j, ct_i, ct_j)
        rp_ji, rs_ji, tp_ji, ts_ji = compute_polarization(n_j, n_i, ct_j, ct_i)

        Ds_i = D_mat(rs_ij, rs_ji, ts_ij, ts_ji)
        Dp_i = D_mat(rp_ij, rp_ji, tp_ij, tp_ji)
        P_i = P_mat(wave, n_i, d[i], ct_i)

        Ts = Ts*P_i*Ds_i
        Tp = Tp*P_i*Dp_i

    r_s = Ts[1,0] / Ts[0,0]
    r_p = Tp[1,0] / Tp[0,0]

    t_s = 1.0 / Ts[0,0]
    t_p = 1.0 / Tp[0,0]

    R_s = modulus_sqrd(r_s)
    R_p = modulus_sqrd(r_p)

    x = n[LAYERS-1] * cos_theta[LAYERS-1]
    y = conjugate(n[LAYERS-1]) * cos_theta[LAYERS-1]
    z = n[0].real * cos_theta[0].real

    T_s = (modulus_sqrd(t_s) * x.real) / z
    T_p = (modulus_sqrd(t_p) * y.real) / z

    R = (R_s + R_p) * 0.5
    T = (T_s + T_p) * 0.5

    return R, T

def sample_interfaces(wave, cos_theta_0, sin_theta_0, n_ior):

    n = [0] * 3
    cos_theta = [0] * 3

    n[0] = sample_IOR(wave, n_ior[0])
    n[1] = sample_IOR(wave, n_ior[1])
    n[2] = sample_IOR(wave, n_ior[2])

    cos_theta[0] = cos_theta_0
    cos_theta[1] = cos_theta_i(n[1], n[0], sin_theta_0)
    cos_theta[2] = cos_theta_i(n[2], n[0], sin_theta_0)

    return n, cos_theta

def imgay(r, film_0_d, film_0_n, film_0_k, base_n, base_k, n, wo):

    d = [INFINITE, film_0_d, INFINITE]

    n_ior = [0] * LAYERS
    n_ior[0] = IOR(np.array((1.0,1.0,1.0)), np.array((0.0,0.0,0.0)))
    n_ior[1] = IOR(film_0_n, film_0_k)
    n_ior[2] = IOR(base_n, base_k)

    Normal = n

    wi = sample_wi(r, wo, Normal)
    h = normalize(wo + wi)

    cos_theta_0 = complex(max(np.dot(h, wo), 1e-6), 0.0)
    sin_theta_0 = sqrt(1.0 - pow(cos_theta_0.real, 2.0))

    lambda_samples = gen_lambda_samples()

    spec_R = [0] * LAMBDA_SAMPLES
    spec_T = [0] * LAMBDA_SAMPLES

    for i in range(LAMBDA_SAMPLES):
        n, cos_theta = sample_interfaces(lambda_samples[i], cos_theta_0, sin_theta_0, n_ior)
        spec_R[i], spec_T[i] = thin_film(lambda_samples[i], n, d, cos_theta)

    #rgb_R = SPEC_to_RGB(spec_R, lambda_samples)
    #rgb_T = SPEC_to_RGB(spec_T, lambda_samples)

    #R = rgb_R
    #T = rgb_T
    return spec_R, spec_T

r = 0.0

film_0_d = 500.0
film_0_n = np.array((2.5,2.5,2.5))
film_0_k = np.array((0.0,0.0,0.0))

base_n = np.array((1.5,1.5,1.5))
base_k =np.array((0.0,0.0,0.0))

n = np.array((0.0,-1.0,0.0))
wo = normalize(np.array((0.0,1.0,0.0)))

plt.title("")

angle_samples = 8

min_angle = 0.0
max_angle = 90.0

angle_step = max_angle / float(angle_samples)

wavelengths = [0] * LAMBDA_SAMPLES
for i in range(LAMBDA_SAMPLES):
    _lambda = (LAMBDA_MIN + (i * LAMBDA_STEP))
    wavelengths[i] = _lambda

for i in range(angle_samples+1):

    angle = math.radians(i * angle_step + min_angle)

    n = np.dot(rotation_matrix([1,0,0], angle), n)

    spec_R, spec_T = imgay(r, film_0_d, film_0_n, film_0_k, base_n, base_k, n, wo)
    plt.plot(wavelengths, spec_R, label=(i * angle_step + min_angle))
    #plt.plot(wavelengths, spec_T, label=r'$T(\lambda)$')

#plt.plot(CIE_L, spec_T, 'red', label=r'$T(\lambda)$')
#plt.plot(wavelengths, D65_coords, 'black', label=r'$D65(\lambda)$')

plt.xlabel(r'$\lambda/nm$')

plt.legend(loc='upper right')

plt.show()