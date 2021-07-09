"""
Simple simulated annealing
No direction
Temperature loop from start to end temperature
Step size decreases with temperature

"""

class pf_gsa:

  t = None
  count = 0

  def run(fit):
    pf_sa.count = pf_sa.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 

    if(pf.fit['sa_loops_t'] == 0 or pf.fit['sa_loops_i'] == 0):
      return 0

    g.pfdata['stage'] = 'Simulated Annealing ' + str(pf_sa.count)
    g.pfdata['stage_brief'] = 'SA' + str(pf_sa.count)

    # Start - use best
    p = g.pfdata['p']['best']
    potential.update(p)
    rss = pf.get_rss(False)
    p_count = len(p)

    # Store Best Parameters
    p_best = numpy.copy(p)
    rss_best = rss

    rss_plot = []
    rss_plot.append([time.time()-t_start, rss_best])

    step = pf.fit['sa_step']  
    for tn in range(pf.fit['sa_loops_t']):
      pgrad = pf_pgradient.run(p_best, g.pfdata['params_var'], True)
      if(pf.fit['sa_loops_t'] == 1):
        pf_sa.t = pf.fit['sa_temp_start']
        f = 1.0
      else:
        pf_sa.t = pf.fit['sa_temp_start'] - tn *  (pf.fit['sa_temp_start'] - pf.fit['sa_temp_end']) / (pf.fit['sa_loops_t'] - 1)
        f = pf.fit['sa_step_factor']**tn
        g.pfdata['stage'] = 'Simulated Annealing ' +str(pf_sa.count) + ' Loop ' + str(tn + 1)
        rss_start_loop = rss_best
      for n in range(pf.fit['sa_loops_i']):
        loop = True
        while(loop):
          p_new = potential.random(p, f * step, False, pgrad)
          rss_new = pf.get_rss(False)
          if(rss_new is not None):
            loop = False
        if(rss_new is not None):
          if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss-rss_new) / pf_sa.t)):
            p = numpy.copy(p_new)
            rss = rss_new
            if(rss_new < rss_best):
              p_best = numpy.copy(p_new)
              rss_best = rss_new
              rss_plot.append([time.time()-t_start, rss_best])
      if(rss_best < rss_start_loop):
        potential.update(p_best)
        pf.get_rss(True)


    pf.summary_line("G_SIMULATED ANNEALING " + str(pf_sa.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("G_SIMULATED_ANNEALING_" + str(pf_sa.count))
    pf_save.rss_plot("G_SIMULATED_ANNEALING_" + str(pf_sa.count), rss_plot)



  





