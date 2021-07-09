import numpy
import matplotlib.pyplot as plt

# SRIM BOOK
def zbl1(r, qa, qb):
  au = (0.8854*0.529)/(qa**0.23+qb**0.23)
  x = r / au
  phi = 0.1818*numpy.exp(-3.2*x)+0.5099*numpy.exp(-0.9423*x)+0.2802*numpy.exp(-0.4029*x)+0.02817*numpy.exp(-0.2016*x)
  return phi * (qa * qb * (1.0 / r))

# LAMMPS
def zbl2(r, qa, qb):
  a = 0.46850/(qa**0.23+qb**0.23)
  x = r / a
  phi = 0.18175*numpy.exp(-3.1998*x)+0.50986*numpy.exp(-0.94229*x)+0.28022*numpy.exp(-0.40290*x)+0.02817*numpy.exp(-0.20162*x)
  pi = 3.1415927
  e = 1.60218e-19
  e0 = 8.854e-12
  jtoev = 6.242e18
  return jtoev * (1 / (4 * pi * e0)) * ((qa * qb * e**2) / (r)) * phi
  
# MENDELEV
def zbl3(r, qa, qb):
  a = 0.4683766/(qa**(2.0/3.0)+qb**(2.0/3.0))
  x = r / a
  phi = 0.1818*numpy.exp(-3.2*x)+0.5099*numpy.exp(-0.9423*x)+0.2802*numpy.exp(-0.4029*x)+0.02817*numpy.exp(-0.2016*x)
  return ((qa * qb) / (r)) * phi


# quantumatk.com
def zbl4(r, qa, qb):
  au = (0.8854*0.529)/(qa**0.23+qb**0.23)
  A = (qa * qb) / (4*3.14159* 8.854e-12)
  x = r / au
  phi = 0.1818*numpy.exp(-3.2*x)+0.5099*numpy.exp(-0.9423*x)+0.2802*numpy.exp(-0.4029*x)+0.02817*numpy.exp(-0.2016*x)
  return phi * (A / r)


# Mendelev Ackland
def zbl5(r, qa, qb):
  au = (0.8854*0.529)/(qa**0.23+qb**0.23)
  x = r / au
  phi = 0.1818*numpy.exp(-3.2*x)+0.5099*numpy.exp(-0.9423*x)+0.2802*numpy.exp(-0.4029*x)+0.02817*numpy.exp(-0.2016*x)
  return (qa * qb * phi) / r


# Mendelev Ackland
def zbl6(r, qa, qb):
  au = (0.8854*0.529)/(qa**0.23+qb**0.23)
  x = r / au
  phi = 0.1818*numpy.exp(-3.2*x)+0.5099*numpy.exp(-0.9423*x)+0.2802*numpy.exp(-0.4029*x)+0.02817*numpy.exp(-0.2016*x)
  return (9734.236 * phi) / r

#print((1/(4*3.14159*8.854e-12)) * (1.60218e-19)**2 * 6.242e18)


d = numpy.zeros((101,6),)
d[:,0] = numpy.linspace(0.01,6.5,101)
d[:,1] = zbl1(d[:,0], 26.0, 26.0)
d[:,2] = zbl2(d[:,0], 26.0, 26.0)
d[:,3] = zbl3(d[:,0], 26.0, 26.0)
d[:,4] = zbl4(d[:,0], 26.0, 26.0)
d[:,5] = zbl5(d[:,0], 26.0, 26.0)

plt.figure(figsize=(12,8))
#plt.plot(d[:,0], d[:,1])
plt.plot(d[:,0], d[:,5])
plt.ylim(-1.0,10.0)
plt.xlabel('')
plt.ylabel('')
plt.title('')
plt.grid(True)
plt.show()
plt.close('all') 
