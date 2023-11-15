import numpy as np
#import pandas as pd
import math
import matplotlib.pyplot as plt

_v = math.sqrt(20.0)
_m = 90.0

np.random.seed(521)
U1 = np.random.uniform(size = 50000)
U2 = np.random.uniform(size = 50000)
R = np.sqrt(-2 * np.log(U1))
Theta = 2 * np.pi * U2
X = R * np.cos(Theta)
Y = R * np.sin(Theta)

_X = _v * X + _m
_Y = _v * Y + _m

print(np.mean(_X), np.mean(_Y))
print(np.var(_X), np.var(_Y))

fig,(ax1,ax2) = plt.subplots(1,2)
temp = ax1.hist(_X)
ax1.set_title("X")
temp = ax2.hist(_Y)
ax2.set_title("Y")
plt.show()