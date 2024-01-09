import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import random

def rand(mu, var):
    u1 = random.random()
    u2 = random.random()
    r = math.sqrt(-2.0 * math.log(u1))
    t = math.pi * 2.0 * u2

    x = r * math.cos(t)

    return math.sqrt(var) * x + mu
'''
l0_sd = 16.0 # lamellae standard deviation of thickness
l0_m = 100.0 # lamellae mean thickness

l1_sd = 16.0 # lamellae standard deviation of thickness
l1_m = 90.0 # lamellae mean thickness

v0 = l0_sd**2 / l0_m**2
v1 = l1_sd**2 / l1_m**2
m = 0.0

L = 1.0#(4.0 * math.pi * 1.56 * l0_m * 1.0) / 600.0

size = 10000

X = [0] * size
Y = [0] * size

for i in range(size):
    X[i] = cmath.exp(complex(0.0, rand(m, v0) * L))
    #X[i] = rand(_m, _v_0)**4
    Y[i] = rand(m, v1)**2

#print(np.mean(np.concatenate((X,Y))))
g = np.mean(X)
print(np.sqrt(g.real*g.real + g.imag*g.imag), np.exp(-L**2 * v0 / 2.0))
'''
'''
print(
    np.mean(np.array([141, 138, 101, 101, 105, 101, 94, 104])), 
    np.mean(np.array([84, 72, 71, 72, 79, 80, 93, 51]))
)
print(
    math.sqrt(np.var(np.array([141, 138, 101, 101, 105, 101, 94, 104]))), 
    math.sqrt(np.var(np.array([84, 72, 71, 72, 79, 80, 93, 51])))
)

fig,(ax1,ax2) = plt.subplots(1,2)
temp = ax1.hist(X)
ax1.set_title("X")
temp = ax2.hist(Y)
ax2.set_title("Y")
#plt.show()
'''