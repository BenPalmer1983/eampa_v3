from f_interp import interp
import numpy


interp.speed_test()






"""
def f(x):
  return 0.03 * x**4 - 1.63 * x**3 + 1.22 * x**2 - 0.0001 * x + 3.7
  
def g(x):
  return 2.3 - 1.2 * x + 0.2 * x**2
  
  
x = numpy.linspace(0,10, 101)
y = f(x[:])


print(x)
print(y)

xi = 2.773
yi = interp.interpolate(xi, x, y)

print(xi, yi, f(xi))



x = numpy.linspace(0,10, 11)
y = f(x[:])

print(x)
print(y)

xi = 0
yi = interp.trap(xi, x, y)
print(xi, yi)

xi = 5
yi = interp.trap(xi, x, y)
print(xi, yi)

xi = 10
yi = interp.trap(xi, x, y)
print(xi, yi)

xi = 5.5
yi = interp.trap(xi, x, y)
print(xi, yi)


x = numpy.linspace(0,10,11)
y = f(x[:])
fill = interp.fill(x, y, 21,4)


print(x)
print(y)
print(fill)
"""


