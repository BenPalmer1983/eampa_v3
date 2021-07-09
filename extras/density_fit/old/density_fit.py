import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import sys
from f_fnc import fnc


class density_fit:

  @staticmethod
  def rss(d, p, pf, f):
    return sum((f(d[:,0], p, pf) - d[:,1])**2)

  @staticmethod
  def dens_a(x, p, pf):
    # f(x) = (x - p[0])**2 
    return (x - p[0])**2 

  @staticmethod
  def dens_b(x, p, pf):
    # Slater 4s
    # f(x) =  (p(1) * r**3 * exp(-1.0D0 * p(2) * r))**2
    return p[0] * x**3 * numpy.exp(-1.0 * p[1] * x)**2 

  @staticmethod
  def dens_c(x, p, pf):
    if(type(x) == numpy.float64):
      return fnc.double_slater_4s(x, p, pf)
    if(type(x) == numpy.ndarray):
      return fnc.double_slater_4s_v(x, p, pf)

  @staticmethod
  def dens_d(x, p, pf):
    if(type(x) == numpy.float64):
      return fnc.slater_4s_cutoff(x, p, pf)
    if(type(x) == numpy.ndarray):
      return fnc.slater_4s_cutoff_v(x, p, pf)
       

  @staticmethod
  def dens_e(x, p, pf):
    if(type(x) == numpy.float64):
      return fnc.double_slater_4s_cutoff(x, p, pf)
    if(type(x) == numpy.ndarray):
      return fnc.double_slater_4s_cutoff_v(x, p, pf)



  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = density_fit.one_space(row)
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
  def fit(file_name, f_type):

    if(f_type == 1):
      f = density_fit.dens_a
      f_size = 1
      f_name = "quadratic_density"
      pf = numpy.asarray([0])
    elif(f_type == 2):
      f = density_fit.dens_b
      f_size = 2
      f_name = "slater_4s"
      pf = numpy.asarray([0])
    elif(f_type == 3):
      f = density_fit.dens_c
      f_size = 4
      f_name = "double_slater_4s"
      pf = numpy.asarray([0])
    elif(f_type == 4):
      f = density_fit.dens_d
      f_size = 2
      f_name = "slater_4s_cutoff"
      pf = numpy.asarray([6.5])
    elif(f_type == 5):
      f = density_fit.dens_e
      f_size = 4
      f_name = "double_slater_4s_cutoff"
      pf = numpy.asarray([6.5])



    pot_file = file_name.split(".")
    pot_file = pot_file[0] + "_a.pot" 
    pot_file_fit = pot_file[0] + "_f.pot" 



    d = density_fit.read_data(file_name)

    p = numpy.zeros((f_size,),)
    p_best = p
    rss_best = None
    for n in range(100):
      for m in range(100):
        p = p_best + density_fit.random_p(f_size, 50 * 0.9999**n)
        rss = density_fit.rss(d, p, pf, f)
        if(rss_best == None or rss < rss_best):
          rss_best = rss
          p_best = p
      #print(rss_best)
        
        
    mf = 0.99
    for n in range(200):
      for m in range(200):
        r = mf**n * (0.5 - numpy.random.rand(len(p)))
        p = p_best + r
        rss = density_fit.rss(d, p, pf, f)
        if(rss_best == None or rss < rss_best):
          rss_best = rss
          p_best = p
    print(rss_best)

    # LMA
    p_best = lma.fit(f, p_best, pf, d[:,0] , d[:,1])
    rss_best = density_fit.rss(d, p_best, pf, f)
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
    plt.show()
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
    fh.close()   







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
  print("python3 density_fit.py densityfile function")
  print()
  print("Function types:")
  print("function 1 = quadratic_density")
  print("function 2 = slater_4s")
  print("function 3 = double_slater_4s")
  print("function 4 = slater_4s_cutoff")
  print("function 5 = double_slater_4s_cutoff")
  print()
  print("Example:")
  print("python3 density_fit.py al_dens.pot 4")
  exit()

filename = str(sys.argv[1]).strip()
func = int(sys.argv[2])
density_fit.fit(filename, func)









