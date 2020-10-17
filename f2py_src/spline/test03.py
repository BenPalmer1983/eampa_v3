import random
import numpy as np
import time
from f_spline import spline
import matplotlib.pyplot as plt


def f(x):
  return 0.01*x**3 + 0.25 * x**2 - 3 * x - 3
def fp(x):
  return 0.03*x**2 + 0.5 * x - 3
def fpp(x):
  return 0.06*x + 0.5

x = np.linspace(0,10,101)
y = f(x)

print(x, y)




arr = spline.spline_array(3, 6, x[:], y[:], 2.0, 8.0, 101, 3)

plt.plot(x[:], y[:])
plt.plot(arr[:,0],arr[:,1])
plt.show()











a = np.zeros((3,),)
b = np.zeros((3,),)


a[0] = 1.0
a[1] = 0.0
a[2] = 1.0

b[0] = 4.0
b[1] = 0.0
b[2] = -1.0


c = spline.spline_ab(1, a, b)

#print(c)


