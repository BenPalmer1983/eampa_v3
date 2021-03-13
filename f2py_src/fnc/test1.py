from f_fnc import fnc
import numpy
import matplotlib.pyplot as plt



p = numpy.zeros((3,),)
p[0] = 0.5
p[1] = 1.2
p[2] = 0.3
x = numpy.linspace(0.0,6.5, 100)
y = fnc.morse_v(x, p)


x = numpy.linspace(-10,10, 100)
y = fnc.heaviside_v(x)



x = numpy.linspace(0,6.5, 100)
p = numpy.zeros((4,),)
p[0] = 2.5
p[1] = 2
p[2] = 0
p[3] = 0
p_fix = numpy.zeros((8,),)
p_fix[4] = 6.5
y = fnc.spline_one_node_v(x, p, p_fix)



x = numpy.linspace(0,6.5, 100)
p = numpy.zeros((12,),)
p[0] = 0
p[1] = 0
p[2] = 0
p[3] = 0
p[4] = 2.5
p[5] = 2
p[6] = 0
p[7] = 0
p[8] = 6.5
p[9] = 0
p[10] = 0
p[11] = 0
y = fnc.spline_n_node_v(x, p)

plt.plot(x, y)
plt.show()


