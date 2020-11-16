class gd:

  precision = 1.0e-6
  max_iterations = 10
  h = 1.0e-8
  momentum = 0.1
  p_in = None
  rss_in = None
  p_out = None
  rss_out = None
  fix = None

  # Central Difference
  def df(p):
    d = numpy.zeros((len(p),),)
    for i in range(len(p)):
      if(gd.fix[i] == 0.0):
        d[i] = 0.0
      else:
        p_f = numpy.copy(p)
        p_f[i] = p_f[i] + gd.h
        p_b = numpy.copy(p)
        p_b[i] = p_b[i] - gd.h
        d[i] = (gd.rss(p_f) - gd.rss(p_b)) / (2 * gd.h)
    return d

  def line_search(p, dp):
    # Back track
    best_rss = gd.rss_in
    best_gamma = 0.0
    gamma = 20.0 * gd.last_gamma 
    if(gamma <1.0e-8):
      gamma <1.0e-6
    loop = True
    n = 0
    while(gamma > 1.0e-12):
      n = n + 1
      p_test = p - gamma * dp
      rss = gd.rss(p_test)
      gamma = 0.5 * gamma
      if(best_rss is None or rss < best_rss):
        p_best = p_test
        best_rss = rss
        best_gamma = gamma
    gd.last_gamma = best_gamma
    return best_rss, best_gamma * dp
    

  # Gradient Descent
  def opt(f_rss, p0, fix=None):
    gd.rss = f_rss   
    gd.last_gamma = 1.0
    p = p0
    
    if(type(fix) == numpy.ndarray):
      if(len(p) == len(fix)):
        for i in range(len(fix)):
          if(fix[i] != 0.0):
            fix[i] = 1.0
      else:
        fix = numpy.zeros((len(p0),),)
        fix[:] = 1.0
    else:
      fix = numpy.zeros((len(p0),),)
      fix[:] = 1.0
    gd.fix = fix  
    gd.p_in = numpy.copy(p)
    gd.rss_in = gd.rss(p)
    
    best_p = numpy.copy(gd.p_in)
    best_rss = gd.rss_in

    n = 0
    last_dp = 0.0
    while(n < gd.max_iterations):
      df = gd.df(p)
      rss, dp = gd.line_search(p, df)
      if(rss > best_rss):
        n = gd.max_iterations
      else:
        p = p - (gd.momentum * last_dp + dp)
        last_dp = dp
        rss = gd.rss(p)
        best_p = numpy.copy(p)
        best_rss = rss
        n = n + 1
        
    gd.p_out = numpy.copy(best_p)
    gd.rss_out = best_rss
    
    return p