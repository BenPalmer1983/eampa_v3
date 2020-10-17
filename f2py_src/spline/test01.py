import random
import numpy as np
import time
from f_spline import spline
import matplotlib.pyplot as plt


def f(x):
  return 0.4 - 1.2 * x + 0.03 * x**3 - 0.002 *x**4

def fp(x):
  return -1.2 + 0.09 * x**2 - 0.008 *x**3
  
def fpp(x):
  return 0.18 * x - 0.024 *x**2
  
def fppp(x):
  return 0.18 - 0.048 * x

n = 201
arr = np.zeros((n,5),)
arr[:,0] = np.linspace(1.0,12.0,n)
arr[:,1] = f(arr[:,0])
arr[:,2] = fp(arr[:,0])
arr[:,3] = fpp(arr[:,0])
arr[:,4] = fppp(arr[:,0])
print(spline.spline_check(arr))
plt.plot(arr[:,0], arr[:,1], marker='o', linestyle='none')
#arr = spline.spline_it(arr, 9)
#arr = spline.spline_make(arr, 9, 1001, 4, 2)

for n in range(6):
  arr = spline.spline_vary_a(arr, 9, 1.2, 0, 2)
  print(spline.spline_check(arr))



plt.plot(arr[:,0], arr[:,1])
plt.show()



#plt.show()

