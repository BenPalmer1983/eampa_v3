import numpy

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter




def plot(data, plot_name, ylim=None):
  plt.figure(figsize=(12,8)) 
  plt.plot(data[:,0], data[:,1])
  plt.xlabel('')
  plt.ylabel('')
  if(ylim != None):
    plt.ylim(ymax = ylim, ymin = -5)
  plt.title('')
  plt.grid(True)
  plt.savefig(plot_name + '.eps', format='eps')
  plt.close('all') 


def H(x):
  if(x<0):
    return 0.0
  return 1.0



def ack_emb(x, p):
  return p[0] * numpy.sqrt(x) + p[1] * x**2 + p[2] * x**4


def density(x, p):
  y = 0.0
  for pn in p:
    y = y + pn[1] * (pn[0] - x)**3 * H(pn[0] - x)
  return y


def pair(r, q, rzbl, b, rb, p):
  if(r == 0.0):
    return 1.0e9

  r1 = rzbl[1]
  r2 = rb[1]

  #rs
  rs = 0.4683766 / (q[0]**(2.0/3.0)+q[1]**(2.0/3.0))

  # ZBL
  y = H(r1 - r) * ((q[0] * q[1]) / r) * e(r / rs)

  # SPLINE
  y = y + H(r2 - r) * H(r - r1) * numpy.exp(b[0] + b[1] * r + b[2] * r**2 + b[3] * r**3)

  for pn in p:
    rk = pn[1]
    y = y + pn[2] * (r - rk)**3 * H(rk - r) * H(r - r2)
    #y = y + pn[2] * (r - rk)**3 * H(rk - r) * H(r - r1)

  #print(rs)
  #print(((q[0] * q[1]) / r) * e(r / rs))
  return y






def pair_test(r, q, rzbl, b, rb, p):
  y = 0.0
  r1 = rzbl[1]
  r2 = rb[1]

  if(r == 0.0):
    y = 1.0e9
  elif(r <= r1):
    rs = 0.4683766 / (q[0]**(2.0/3.0)+q[1]**(2.0/3.0))
    y = H(r1 - r) * ((q[0] * q[1]) / r) * e(r / rs)
  elif(r > r1 and r < r2):
    rs = 0.4683766 / (q[0]**(2.0/3.0)+q[1]**(2.0/3.0))
    
    xa = r1
    ya = H(r1 - xa) * ((q[0] * q[1]) / xa) * e(xa / rs)
    xb = r2
    yb = 0.0

    for pn in p:
      rk = pn[1]
      yb = yb + pn[2] * (xb - rk)**3 * H(rk - xb) * H(xb - r2)
  
    xmat = numpy.zeros((4,4,),)
    ymat = numpy.zeros((4,),) 
    
    xmat[0,0] = 1
    xmat[0,1] = xa
    xmat[0,2] = xa**2
    xmat[0,3] = xa**3
    xmat[1,0] = 0
    xmat[1,1] = 1
    xmat[1,2] = 2 * xa
    xmat[1,3] = 3 * xa**2
    xmat[2,0] = 1
    xmat[2,1] = xb
    xmat[2,2] = xb**2
    xmat[2,3] = xb**3
    xmat[3,0] = 0
    xmat[3,1] = 1
    xmat[3,2] = 2 * xb
    xmat[3,3] = 3 * xb**2
    
    ymat[0] = ya
    ymat[1] = 0.0
    ymat[2] = yb
    ymat[3] = 0.0
    
    c = numpy.linalg.solve(xmat, ymat)
    y = c[0] + c[1] * r + c[2] * r**2 + c[3] * r**3
    
  elif(r >= r2):
    for pn in p:
      rk = pn[1]
      y = y + pn[2] * (r - rk)**3 * H(rk - r) * H(r - r2)

  return y


def e(x):
  return 0.1818 * numpy.exp(-3.2 * x) + 0.5099 * numpy.exp(-0.9423 * x) + 0.2802 * numpy.exp(-0.4029 * x) + 0.02817 * numpy.exp(-0.2016 * x)


def make(pot, pot_size=1001, pmin=0.0, pmax = 6.5, rhomin=0.0, rhomax = 6.5, fmin=0.0, fmax=1.0):

  # Pair
  v = numpy.zeros((pot_size, 2),)
  v[:,0] = numpy.linspace(rhomin, rhomax, pot_size)
  for n in range(pot_size):
    v[n,1] = pair_test(v[n,0], pot['Q'], pot['RZBL'], pot['B'], pot['RSPLINE'], pot['PAIR'])

  # Density
  rho = numpy.zeros((pot_size, 2),)
  rho[:,0] = numpy.linspace(rhomin, rhomax, pot_size)
  for n in range(pot_size):
    rho[n,1] = density(rho[n,0], pot['RHO'])

  # Embedding
  f = numpy.zeros((pot_size, 2),)
  f[:,0] = numpy.linspace(fmin, fmax, pot_size)
  f[:,1] = ack_emb(f[:,0], pot['F'])

  # Plot
  plot(v, 'pair', 500)
  plot(rho, 'dens')
  plot(f, 'embe')
  
  fh = open("fe_pair.pot", 'w')
  for n in range(len(v)):
    fh.write(str(v[n,0]) + "   " + str(v[n,1]) + "\n")
  fh.close()
  
  
  fh = open("fe_dens.pot", 'w')
  for n in range(len(rho)):
    fh.write(str(rho[n,0]) + "   " + str(rho[n,1]) + "\n")
  fh.close()
  
  
  fh = open("fe_embe.pot", 'w')
  for n in range(len(f)):
    fh.write(str(f[n,0]) + "   " + str(f[n,1]) + "\n")
  fh.close()
  
  
  
  return v, rho, f



def load(pname):
  pot = {
         'Q': [],
         'B': [],
         'PAIR': [],
         'RHO': [],
         'F': [],
         'RZBL': [],
         'RSPLINE': [],
         'RHO': [],
         'F': [],
         }


  fh = open(pname, 'r')
  for row in fh:
    row = row.strip().split(' ')

    if(row[0][0] == 'Q'):
      pot['Q'].append(float(row[2]))
    elif(row[0][0:5] == 'RZBL'):
      pot['RZBL'].append(float(row[2]))
    elif(row[0][0] == 'B'):
      pot['B'].append(float(row[2]))
    elif(row[0][0:11] == 'RZBLSPLINE'):
      pot['RSPLINE'].append(float(row[2]))
    elif(row[0][0:11] == 'P'):
      pot['PAIR'].append([float(row[2]),float(row[3]),float(row[4])])
    elif(row[0][0:4] == 'RHO'):
      pot['RHO'].append([float(row[2]),float(row[3])])
    elif(row[0][0] == 'F'):
      pot['F'].append(float(row[2]))
  

  fh.close()
  return pot


numpy.set_printoptions(threshold=10000)

pot = load("parameters.in")

make(pot)


