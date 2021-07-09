class pf_pgradient:

  def run(p, p_var=None, direction=False):
    if(type(p_var) == type(None)):
      p_var = numpy.zeros((len(p),),)
      p_var[:] = 1.0

    potential.update(p)
    rss = pf.get_rss(False, True)
    grad = numpy.zeros((len(p),), dtype=numpy.float64,)

    for i in range(len(p)):
      if(p_var[i] == 0):
        grad[i] = 0.0
      else:
        # Forward
        p_f = numpy.copy(p)
        h = 1.0e-8
        p_f[i] = p_f[i] + h
        potential.update(p_f)
        rss_f = pf.get_rss(False, True)
        # Backward
        p_b = numpy.copy(p)
        h = 1.0e-8
        p_b[i] = p_b[i] - h
        potential.update(p_b)
        rss_b = pf.get_rss(False, True)
      
        grad[i] = (rss_f - rss_b) / (2.0 * h)
    if(direction == True):
      for i in range(len(grad)):
        if(grad[i] != 0):
          grad[i] = grad[i] / abs(grad[i])
    return grad
