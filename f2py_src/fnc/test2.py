from f_fnc import fnc
import numpy



    
x = numpy.linspace(0.5, 6.5, 101) 
p = numpy.asarray([0.1,-0.1,0.001,0.0001])  
p_fixed = numpy.asarray([26.0,26.0,1.0,2.0,3.0,4.0,5.0,6.5])


x = 0.5
y = fnc.pair_spline(x, p, p_fixed)
print(x, y)
x = 1.5
y = fnc.pair_spline(x, p, p_fixed)
print(x, y)
x = 10.0
y = fnc.pair_spline(x, p, p_fixed)
print(x, y)
