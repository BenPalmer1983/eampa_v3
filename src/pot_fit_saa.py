"""
Simple simulated annealing
No direction
Temperature loop from start to end temperature
Step size decreases with temperature

"""

class pf_saa:

  t = None
  count = 0

  def run():
    pf_saa.count = pf_saa.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 

    #if(pf.fit['sa_loops_t'] == 0 or pf.fit['sa_loops_i'] == 0):
    #  return 0

    # Start - use best
    p = g.pfdata['p']['best']
    pf_potential.update(p)
    rss = pf.get_rss(False)
    p_count = len(p)

    # Store Best Parameters
    p_best = numpy.copy(p)
    rss_best = rss
    g.pfdata['stage'] = 'Simulated Annealing Advanced ' + str(pf_saa.count) + ' Line Search'

    f = 100.0
    f_best = f
    p_grad = pf_pgradient.run(p_best)
    for n in range(100):
      p_new = numpy.copy(p_best)
      p_new = p_new - f * p_grad
      pf_potential.update(p_new)
      rss_new = pf.get_rss(False)
      if(rss_new < rss_best):
        rss_best = rss_new 
        f_best = f
      f = 0.5 * f
 
    
    g.pfdata['stage'] = 'Simulated Annealing Advanced ' + str(pf_saa.count) + ' Anneal'

    for i in range(10):
      for n in range(50):
        p_new = numpy.copy(p_best)
        p_new = p_new - 0.1 * f_best * p_grad * numpy.random.rand(p_count)   

        pf_potential.update(p_new)
        rss_new = pf.get_rss(False)

        if(rss_new < rss_best):
          p_best = numpy.copy(p_new)
          rss_best = rss_new
      p_grad = pf_pgradient.run(p_best)
   

    exit()

    """
    step = pf.fit['sa_step']  
    for tn in range(pf.fit['sa_loops_t']):
      if(pf.fit['sa_loops_t'] == 1):
        pf_sa.t = pf.fit['sa_temp_start']
        f = 1.0
      else:
        pf_sa.t = pf.fit['sa_temp_start'] - tn *  (pf.fit['sa_temp_start'] - pf.fit['sa_temp_end']) / (pf.fit['sa_loops_t'] - 1)
        f = pf.fit['sa_step_factor']**tn
        g.pfdata['stage'] = 'Simulated Annealing Loop ' + str(tn + 1)
        rss_start_loop = rss_best
      for n in range(pf.fit['sa_loops_i']):
        p_new = numpy.copy(p)
        p_new = p_new + f * step * (0.5 - numpy.random.rand(p_count))
        pf_potential.update(p_new)
        rss_new = pf.get_rss(False)

        if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss-rss_new) / pf_sa.t)):
          p = numpy.copy(p_new)
          rss = rss_new
          if(rss_new < rss_best):
            p_best = numpy.copy(p_new)
            rss_best = rss_new
      if(rss_best < rss_start_loop):
        pf_potential.update(p_best)
        rss_new = pf.get_rss(True)
    """

    pf.summary_line("SIMULATED ANNEALING " + str(pf_random.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])



  





