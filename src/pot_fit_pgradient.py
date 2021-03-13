class pf_pgradient:

  def run(p):
    pf_potential.update(p)
    rss = pf.get_rss(False)
    grad = numpy.zeros((len(p),), dtype=numpy.float64,)

    for i in range(len(p)):
      # Forward
      p_f = numpy.copy(p)
      h = 1.0e-10 * p_f[i]
      if(p_f[i] == 0.0):
        h = 1.0e-10
      p_f[i] = p_f[i] + h
      pf_potential.update(p_f)
      rss_f = pf.get_rss(False)
      # Backward
      p_b = numpy.copy(p)
      h = 1.0e-10 * p_b[i]
      if(p_b[i] == 0.0):
        h = 1.0e-10
      p_b[i] = p_b[i] - h
      pf_potential.update(p_b)
      rss_b = pf.get_rss(False)
      
      grad[i] = (rss_f - rss_b) / (2.0 * h)

    return grad
