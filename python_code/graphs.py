def f_thickness(wavelength, thickness, film_a, film_b):

    R = [0] * len(wavelength)

    for i in range(len(wavelength)):

        film_a.d = thickness[i]

        l0_sd = 16.0 # lamellae standard deviation of thickness
        film_a.v = l0_sd**2 / film_a.d**2

        sin_theta = sin(radians(15.0))
        wave = wavelength[i]

        phase_a = phase(film_a, wave, sin_theta)
        phase_b = phase(film_b, wave, sin_theta)

        alpha_a = cgf(film_a, phase_a)
        alpha_b = cgf(film_b, phase_b)

        phi = np.exp(-2.0*(alpha_a + alpha_b))
        psi = np.exp(-1.0*(alpha_a + alpha_b))

        x = (np.exp(-alpha_a) + np.exp(-alpha_b))/(1.0 + psi)
        y = (np.exp(-alpha_a) - np.exp(-alpha_b))/(1.0 - psi)

        M = (phase_a + phase_b) * .5
        N = (phase_a - phase_b) * .5

        numer = 1 - cos(M)*cos(N)*x + sin(M)*sin(N)*y
        denom = (1 + phi - 2.0 * cos(2.0*M) * psi)

        R[i] = (1.0 - phi) * numer / denom

    return R

def f_sd(wavelength, sd, film_a, film_b):

    R = [0] * len(wavelength)

    for i in range(len(wavelength)):

        film_a.v = sd[i]**2 / film_a.d**2

        sin_theta = sin(radians(15.))
        wave = wavelength[i]

        phase_a = phase(film_a, wave, sin_theta)
        phase_b = phase(film_b, wave, sin_theta)

        alpha_a = cgf(film_a, phase_a)
        alpha_b = cgf(film_b, phase_b)

        phi = np.exp(-2.0*(alpha_a + alpha_b))
        psi = np.exp(-1.0*(alpha_a + alpha_b))

        x = (np.exp(-alpha_a) + np.exp(-alpha_b))/(1.0 + psi)
        y = (np.exp(-alpha_a) - np.exp(-alpha_b))/(1.0 - psi)

        M = (phase_a + phase_b) * .5
        N = (phase_a - phase_b) * .5

        numer = 1 - cos(M)*cos(N)*x + sin(M)*sin(N)*y
        denom = (1 + phi - 2.0 * cos(2.0*M) * psi)

        maximum = 1.0#(1.0+psi)/(1.0-psi)

        R[i] = (1.0 - phi) * numer / (denom*maximum)

    return R

def f_angle(wavelength, angle, film_a, film_b):

    R = [0] * len(wavelength)

    for i in range(len(wavelength)):

        sin_theta = sin(radians(angle[i]))
        wave = wavelength[i]

        phase_a = phase(film_a, wave, sin_theta)
        phase_b = phase(film_b, wave, sin_theta)

        alpha_a = cgf(film_a, phase_a)
        alpha_b = cgf(film_b, phase_b)

        phi = np.exp(-2.0*(alpha_a + alpha_b))
        psi = np.exp(-1.0*(alpha_a + alpha_b))

        x = (np.exp(-alpha_a) + np.exp(-alpha_b))/(1.0 + psi)
        y = (np.exp(-alpha_a) - np.exp(-alpha_b))/(1.0 - psi)

        M = (phase_a + phase_b) * .5
        N = (phase_a - phase_b) * .5

        numer = 1 - cos(M)*cos(N)*x + sin(M)*sin(N)*y
        denom = (1 + phi - 2.0 * cos(2.0*M) * psi)

        R[i] = (1.0 - phi) * numer / denom

    return R

# ------------------------------- Red Schiller ------------------------------- #

l0_m = 176.6 # lamellae mean thickness
l0_sd = 16.0 # lamellae standard deviation of thickness
l0_v = l0_sd**2 / l0_m**2
film_a = Film(l0_m, l0_v, 1.560, 0.0)

l1_m = 100.4 # lamellae mean thickness
l1_sd = 16.0 # lamellae standard deviation of thickness
l1_v = l1_sd**2 / l1_m**2
film_b = Film(l1_m, l1_v, 1.56, 0.0)

# ------------------------------- Blue Schiller ------------------------------ #
'''
l0_m = 72.5 # lamellae mean thickness
l0_sd = 16.0 # lamellae standard deviation of thickness
l0_v = l0_sd**2 / l0_m**2
film_a = Film(l0_m, l0_v, 1.560, 0.0)

l1_m = 65.1 # lamellae mean thickness
l1_sd = 16.0 # lamellae standard deviation of thickness
l1_v = l1_sd**2 / l1_m**2
film_b = Film(l1_m, l1_v, 1.560, 0.0)
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

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
max_sd = 60.0
sd_step = float(max_sd - min_sd) / float(CIE_SAMPLES-1)
sd = [(x * sd_step + min_sd) for x in range(CIE_SAMPLES)]

X, Y = np.meshgrid(wavelength, sd)
zs = np.array(f_sd(np.ravel(X), np.ravel(Y), film_a, film_b))

Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel(r'$\lambda/nm$')
ax.set_ylabel(r'$\theta$')
ax.set_zlabel('Z Label')

plt.xlim([CIE_MIN, CIE_MAX])

#ax.set_zlim([0, 1])

plt.gca().invert_yaxis()

plt.show()