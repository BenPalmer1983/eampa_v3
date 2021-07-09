import sys
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from eampa_lib.f_fnc import fnc


class embe_fit:



  @staticmethod
  def rss(d, p, pf, f):
    y = fnc.fv(embe_fit.fn, d[:,0], p, pf)
    return sum((y - d[:,1])**2)


  @staticmethod
  def f(x, p, pf):
    return fnc.f(embe_fit.fn, x, p, pf)
    
  @staticmethod
  def fv(x, p, pf):
    return fnc.fv(embe_fit.fn, x, p, pf)


  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = embe_fit.one_space(row)
        f = row.split(" ")
        
        line = []
        for fn in f:
          line.append(float(fn))
        
        d.append(line)

    d = numpy.asarray(d)
    return d
          

  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out


  @staticmethod
  def random_p(f_size, f=1.0):
    return (f * (0.5-numpy.random.rand(f_size)))

  @staticmethod
  def fit(file_name, dirout):
  
    embe_fit.make_dir(dirout)
  

    numpy.random.seed(1)
  
    embe_fit.fn = 'embedding_g'
    pf = numpy.asarray([0.1])
    p = numpy.zeros((2,),)
    p[0] = -7.8
    p[1] = -1.5
    
    
    f = embe_fit.f
    fv = embe_fit.fv
    

    pot_file = file_name.split(".")
    pot_file = dirout + '/' + pot_file[0] + "_a.pot" 
    pot_file_fit = dirout + '/' + pot_file[0] + "_f.pot" 
    plot_file_fit = dirout + '/' + pot_file[0] + "_plot.eps" 



    d = embe_fit.read_data(file_name)
    
    p_best = numpy.copy(p)
    rss_best = embe_fit.rss(d, p, pf, f)

        
    d_fit = numpy.copy(d)
    d_fit[:,1] = f(d_fit[:,0], p_best, pf)
       
    for n in range(10000):
      p = numpy.copy(p_best) 
      p[:] = 100.0 * (0.5-numpy.random.rand(len(p)))
      rss = embe_fit.rss(d, p, pf, f)

      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p

    
    for n in range(100000):
      p = numpy.copy(p_best)
      p[:] = p[:] + embe_fit.random_p(len(p), 0.9999**n)
      rss = embe_fit.rss(d, p, pf, f)
      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p
        print(rss_best)
    print(rss_best)
    print(p)
    
    d_fit = numpy.copy(d)
    d_fit[:,1] = fv(d_fit[:,0], p_best, pf)


    plt.figure(figsize=(12,8))
    plt.plot(d[:,0], d[:,1])
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig(plot_file_fit, format='eps')
    plt.show()
    plt.close('all') 
        
        
        
    fh = open(pot_file, "w")
    fh.write("#TYPE " + str(embe_fit.fn) + " \n")
    fh.write("#P ")
    for i in range(len(p_best)):
      fh.write(str(p_best[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(pf)):
      fh.write(str(pf[i]) + " ")
    fh.write("\n")
    fh.write("#VR 0.1\n")
    fh.close()   
        
        
        
        
    exit()
        
    x = numpy.linspace(0.0,20.0)
    y = fv(x, p_best, pf)


    plt.figure(figsize=(12,8))
    plt.plot(d[:,0], d[:,1])
    plt.plot(x, y)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    #plt.savefig('fit2.eps', format='eps')
    plt.show()
    plt.close('all')     
        
    """
    p = numpy.zeros((f_size,),)
    p_best = p
    rss_best = None
    for n in range(10000):
      p = p_best + embe_fit.random_p(f_size, 0.9999**n)
      rss = embe_fit.rss(d, p, pf, f)
      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p
    print(rss_best)
    
    
    
    # LMA
    p_best = lma.fit(f, p_best, pf, d[:,0] , d[:,1])
    rss_best = embe_fit.rss(d, p_best, pf, f)
    print(rss_best)
    
    p_best = lma.fit_ng(f, p_best, pf, d[:,0] , d[:,1])
    rss_best = embe_fit.rss(d, p_best, pf, f)
    print(rss_best)
    
    

    d_fit = numpy.copy(d)
    d_fit[:,1] = f(d_fit[:,0], p_best, pf)


    plt.figure(figsize=(12,8))
    plt.plot(d[:,0], d[:,1])
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('fit.eps', format='eps')
    plt.close('all') 



    fh = open(pot_file_fit, "w")
    for i in range(len(d_fit)):
      fh.write(str(d_fit[i,0]) + " " + str(d_fit[i,0]) + "\n")
    fh.close() 


 
    fh = open('out.pot', "w")
    fh.write("#TYPE " + f_name + " \n")
    fh.write("#P ")
    for i in range(len(p_best)):
      fh.write(str(p_best[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(pf)):
      fh.write(str(pf[i]) + " ")
    fh.write("\n")
    fh.write("#VR 0.1\n")
    fh.close()   
    
    """
    
    
  @staticmethod
  def make_dir(dir):
    dirs = dir.split("/")
    try:
      dir = ''
      for i in range(len(dirs)):
        dir = dir + dirs[i]
        if(not os.path.exists(dir) and dir.strip() != ''):
          os.mkdir(dir) 
        dir = dir + '/'
      return True
    except:
      return False



class lma:

  @staticmethod
  def fit(f, p, pf, x, y, conv=1.0e-12, outer=10, inner=10):
    p_best = numpy.copy(p)
    rss_best = lma.rss(f, p_best, pf, x, y)
    
    loop = True
    cycle = 0
    while(loop):

      R = lma.residual(f, p_best, pf, x, y)        # R
      J = lma.jacobian(f, p_best, pf, x, y)        # J
      JT = numpy.transpose(J)
      JTJ = numpy.matmul(JT, J)
      JTR = numpy.matmul(JT, R)  
      lcut = lma.getlambda(J)
      l = lcut
      rss_last = rss_best
      
      for n in range(inner):
        A = JTJ + l * numpy.diag(JTJ)
        dp = numpy.linalg.solve(A, -JTR)
        p_new = p_best + dp
        rss = lma.rss(f, p_new, pf, x, y)
        
        if(rss < rss_best):
          l = l * 0.2
          if(l < lcut):
            l = lcut
          p_best = numpy.copy(p_new)
          rss_best = rss
        else:
          l = l * 1.5
                 
      
      # Loop counter
      cycle = cycle + 1
      if(cycle > outer or abs(rss_last - rss_best) < conv):
        loop = False
      
    return p_best
    
    
  @staticmethod
  def fit_ng(f, p, pf, x, y):
    p_best = numpy.copy(p)
    rss_best = lma.rss(f, p_best, pf, x, y)
    
    loop = True
    cycle = 0
    while(loop):

      R = lma.residual(f, p, pf, x, y)        # R
      J = lma.jacobian(f, p, pf, x, y)        # J
      JT = numpy.transpose(J)
      JTJ = numpy.matmul(JT, J)
      JTR = numpy.matmul(JT, R)  
 
      dp = numpy.linalg.solve(JTJ, -JTR)
      p = p + dp
      rss = lma.rss(f, p, pf, x, y)
      if(rss < rss_best):
        rss_best = rss
        p_best = numpy.copy(p)
      else:
        loop = False
      cycle = cycle + 1
      if(cycle > 10):
        loop = False
    return p_best
    
    

  @staticmethod
  def residual(f, p, pf, x, y):
    # Calculate residual
    return f(x[:], p, pf) - y[:] 


  @staticmethod
  def jacobian(f, p, pf, x, y):
    dl = len(x)
    pl = len(p)
    J = numpy.zeros((dl, pl), dtype=numpy.float64)    
    h = 0.0001 * ((x[-1] - x[0]) / len(x))

    for i in range(0, dl):
      for j in range(0, pl):
        # Reset parameters
        p_backward = numpy.copy(p)
        p_forward = numpy.copy(p)
        # Vary jth parameter
        p_backward[j] = p_backward[j] - 0.5 * h  
        p_forward[j] = p_forward[j] + 0.5 * h      
        # Calc J matrix
        J[i,j] = (f(x[i], p_forward, pf) - f(x[i], p_backward, pf)) / h
    return J


  @staticmethod
  def rss(f, p, pf, x, y):
    return sum((f(x[:], p, pf) - y[:])**2)


  @staticmethod
  def getlambda(J):
    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_inv = numpy.linalg.inv(JTJ)
    trace = numpy.trace(JTJ_inv)
    l = abs(trace**(-1))
    return l


    
if(len(sys.argv) != 3):
  print("How to use:")
  print("python3 embe_ackland.py embedfile outdir")
  print()
  print()
  print("Example:")
  print("python3 embe_ackland.py Pd_plot.embed pdembed")
  exit()

filename = str(sys.argv[1]).strip()
outdir = str(sys.argv[2]).strip()
embe_fit.fit(filename, outdir)









