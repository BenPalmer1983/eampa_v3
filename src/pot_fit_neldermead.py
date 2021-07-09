


"""
Simple simulated annealing
No direction
Temperature loop from start to end temperature
Step size decreases with temperature

"""

class pf_neldermead:

  count = 0
  p_full = None

  def run(fit):
    for repeat in range(fit['repeat']):
      pf_neldermead.run_nm(fit)

  def get_rss(ps):
    p = pf_neldermead.grow_p(ps)
    potential.update(p)
    rss = pf.get_rss(False)
    return rss

  def shrink_p(p):
    p_out = []
    pdiff = potential.pot['p_upper'][:] - potential.pot['p_lower'][:]
    for i in range(len(pdiff)):
      if(pdiff[i] != 0.0):
        p_out.append(p[i])
    p_out = numpy.asarray(p_out)
    return p_out

  def grow_p(psmall):
    p_out = []
    pdiff = potential.pot['p_upper'][:] - potential.pot['p_lower'][:]
    n = 0
    for i in range(len(pdiff)):
      if(pdiff[i] != 0.0):
        p_out.append(psmall[n])
        n = n + 1
      else:
        p_out.append(pf_neldermead.p_full[i])
    p_out = numpy.asarray(p_out)
    return p_out
    

  def run_nm(fit):
    pf_neldermead.count = pf_neldermead.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 


    g.pfdata['stage'] = 'Nelder-Mead ' + str(pf_neldermead.count)
    g.pfdata['stage_brief'] = 'NM' + str(pf_neldermead.count)
  
    pf_neldermead.p_full = numpy.copy(g.pfdata['p']['best'])
    p = numpy.copy(g.pfdata['p']['best'])

    # Shrink to non-fixed parameters
    ps = pf_neldermead.shrink_p(p)

    res = minimize(pf_neldermead.get_rss, ps, method='nelder-mead',
              options={'xatol':  float(fit['nm_xatol']), 'maxiter': int(fit['nm_maxiter']), })
    ps = numpy.copy(res['x'])

    p = pf_neldermead.grow_p(ps)
    potential.update(p)
    rss = pf.get_rss(False)


    pf.summary_line("NELDER MEAD " + str(pf_neldermead.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("NELDER_MEAD_" + str(pf_neldermead.count))
    #pf_save.rss_plot("NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)
    #pf.save_rss_plot_data(t_start, "NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)







 