import numpy
from eampa_lib.f_fnc import fnc
import matplotlib.pyplot as plt


print("Fixed End")

pf = numpy.zeros((10,),)
p = numpy.zeros((7,),)
x = numpy.linspace(0.0,7.0,1001)


pf[0] = 0.0
pf[1] = 1.0
pf[2] = 2.0
pf[3] = 3.0
pf[4] = 4.0
pf[5] = 5.0
pf[6] = 6.0
pf[7] = 7.0
pf[8] = 7.0
pf[9] = -10.0

p[0] = 0.0
p[1] = 1.0
p[2] = 2.0
p[3] = 1.0
p[4] = 4.0
p[5] = 3.0
p[6] = 6.0



y = fnc.cubic_knot_spline_fixed_end_v(x, p, pf)


plt.plot(x,y)
plt.show()




print("Points")

pf = numpy.zeros((8,),)
p = numpy.zeros((8,),)
x = numpy.linspace(0.0,7.0,1001)


pf[0] = 0.0
pf[1] = 1.0
pf[2] = 2.0
pf[3] = 3.0
pf[4] = 4.0
pf[5] = 5.0
pf[6] = 6.0
pf[7] = 7.0

p[0] = 0.0
p[1] = 1.0
p[2] = 2.0
p[3] = 4.0
p[4] = 2.0
p[5] = 1.0
p[6] = 0.5
p[7] = 0.0



y = fnc.cubic_knot_spline_v(x, p, pf)


plt.plot(x,y)
plt.show()





pf = numpy.zeros((11,),)
p = numpy.zeros((5,),)
x = numpy.linspace(0.0,7.0,1001)


pf[0] = 2.0
pf[1] = 3.0
pf[2] = 4.0
pf[3] = 5.0
pf[4] = 6.0
pf[5] = 13.0    # qa
pf[6] = 13.0    # qb
pf[7] = 1.0     # rzbl
pf[8] = 6.5     # fixed end x
pf[9] = 0.0     # fixed end v(x)
pf[10] = 0.0    # fixed end v'(x) 

p[0] = -0.1
p[1] = 0.1
p[2] = 0.05
p[3] = -0.05
p[4] = 0.01





y = fnc.cubic_knot_spline_fixed_end_pair_v(x, p, pf)



plt.ylim(-1,10)
plt.plot(x,y)
plt.show()





