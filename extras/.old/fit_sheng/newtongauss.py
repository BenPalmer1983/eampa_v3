import numpy


class newtongauss:

  # f   is the vectorised function
  # fs  is scalar version of function
  # p   adjustable parameters of function
  # pf  fixed parameters of function
  # x array
  # y array

  @staticmethod
  def fit(f, p, pf, x, y, fs = None):
    if(fs == None):
      fs = f
      
    p_best = p
    rss_best = newtongauss.rss(f, p, pf, x, y)
    
    loop = True
    cycle = 0
    while(loop):

      R = newtongauss.residual(f, p, pf, x, y)        # R
      J = newtongauss.jacobian(f, fs, p, pf, x, y)        # J

      JT = numpy.transpose(J)
      JTJ = numpy.matmul(JT, J)
      JTR = numpy.matmul(JT, R)  
 
      dp = numpy.linalg.solve(JTJ, -JTR)
      p = p + dp
      if(newtongauss.rss(f, p, pf, x, y) < rss_best):
        rss_best = newtongauss.rss(f, p, pf, x, y)
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
  def jacobian(f, fs, p, pf, x, y):
    dl = len(x)
    pl = len(p)
    J = numpy.zeros((dl, pl), dtype=numpy.float64)  
    h = 0.000001

    for i in range(0, dl):
      for j in range(0, pl):
        # Reset parameters
        p_b = numpy.copy(p)
        p_f = numpy.copy(p)
        # Vary jth parameter
        p_b[j] = p_b[j] - 0.5 * h  
        p_f[j] = p_f[j] + 0.5 * h      
        # Calc J matrix
        
        J[i,j] = (fs(x[i], p_f, pf) - fs(x[i], p_b, pf)) / h

    return J


  @staticmethod
  def rss(f, p, pf, x, y):
    return sum((f(x[:], p, pf) - y[:])**2)




