


"""
Simple simulated annealing
No direction
Temperature loop from start to end temperature
Step size decreases with temperature

"""

class pf_bfgs:

  count = 0
  p_full = None

  def run(fit):
    for repeat in range(fit['repeat']):
      pf_bfgs.run_nm(fit)

  def get_rss(ps):
    p = pf_bfgs.grow_p(ps)
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
        p_out.append(pf_bfgs.p_full[i])
    p_out = numpy.asarray(p_out)
    return p_out
    

  def run_nm(fit):
    pf_bfgs.count = pf_bfgs.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 


    g.pfdata['stage'] = 'BFGS ' + str(pf_bfgs.count)
    g.pfdata['stage_brief'] = 'BFGS' + str(pf_bfgs.count)
  
    pf_bfgs.p_full = numpy.copy(g.pfdata['p']['best'])
    p = numpy.copy(g.pfdata['p']['best'])

    # Shrink to non-fixed parameters
    ps = pf_bfgs.shrink_p(p)

    res = minimize(pf_bfgs.get_rss, ps, method='BFGS',
              options={'gtol':  float(fit['bfgs_gtol']), 'maxiter': int(fit['bfgs_maxiter']), })
    ps = numpy.copy(res['x'])

    p = pf_bfgs.grow_p(ps)
    potential.update(p)
    rss = pf.get_rss(False)


    pf.summary_line("BFGS " + str(pf_bfgs.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("BFGS_" + str(pf_bfgs.count))







 