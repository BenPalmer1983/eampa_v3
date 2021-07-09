from f_neighbourlist import neighbourlist
import numpy

rcut = 6.5
a0 = 4.04

uv = numpy.zeros((3,3,),)
uv[0,0] = 1.0
uv[1,1] = 1.0
uv[2,2] = 1.0

copies = numpy.zeros((3,),)
copies[0] = 4
copies[1] = 4
copies[2] = 4

labels = numpy.zeros((4,),dtype=numpy.int32)

coords = numpy.zeros((4,3,),)
coords[0,0] = 0.0
coords[0,1] = 0.0
coords[0,2] = 0.0
coords[1,0] = 0.5
coords[1,1] = 0.5
coords[1,2] = 0.0
coords[2,0] = 0.0
coords[2,1] = 0.5
coords[2,2] = 0.5
coords[3,0] = 0.5
coords[3,1] = 0.0
coords[3,2] = 0.5


neighbourlist.init()
for i in range(100):
  uv[0,0] = 1.0
  uv[1,1] = 1.0
  uv[2,2] = 1.0 + i * 0.01
  neighbourlist.add_config(rcut, a0, uv, copies, labels, coords)


neighbourlist.build()

print(neighbourlist.cc)
print(neighbourlist.nl_key[0:neighbourlist.cc, 2])
#print(neighbourlist.coords_real[0:3, 0:3])
#print(neighbourlist.nl_key[0:3, 0:3])

#print(neighbourlist.nl_r[0:1000])
