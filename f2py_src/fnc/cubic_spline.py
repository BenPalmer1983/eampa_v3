from f_fnc import fnc
import numpy
import matplotlib.pyplot as plt


p = numpy.asarray([0.01,0.0,-0.01,0.0])
pf = numpy.asarray([2.0,3.0,4.0,5.0,26.0,26.0,0.8,1.8,1])


x = numpy.linspace(0.0, 7.0, 1001)
y = fnc.fv("cubic_spline_zbl", x, p, pf)

plt.plot(x,y,'k')
plt.ylim(-5.0, 10.0)
plt.show()
