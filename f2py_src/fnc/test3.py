from f_fnc import fnc
import numpy

p = numpy.asarray([5.0])
pf = numpy.asarray([0.0])



y = fnc.f("fs_embedding",2.0, p, pf)
print(y)

