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

nodes = spline.get_nodes(x, y, 2.0, 8.5, 7)











