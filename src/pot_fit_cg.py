
# Basin Hopping

class pf_cg:

  count = 0
  p_full = None

  def run(fit):
    for repeat in range(fit['repeat']):
      pf_cg.run_cg(fit)

  def get_rss(ps):
    p = pf_cg.grow_p(ps)
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
        p_out.append(pf_cg.p_full[i])
    p_out = numpy.asarray(p_out)
    return p_out
    

  def run_cg(fit):
    pf_cg.count = pf_cg.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 


    g.pfdata['stage'] = 'Conjugate-Gradient ' + str(pf_cg.count)
    g.pfdata['stage_brief'] = 'CG' + str(pf_cg.count)
  
    pf_cg.p_full = numpy.copy(g.pfdata['p']['best'])
    p = numpy.copy(g.pfdata['p']['best'])

    # Shrink to non-fixed parameters
    ps = pf_cg.shrink_p(p)
    res = minimize(pf_cg.get_rss, ps, method='CG')
    ps = numpy.copy(res['x'])
    p = pf_cg.grow_p(ps)
    potential.update(p)
    rss = pf.get_rss(False)


    pf.summary_line("CONJUGATE GRADIENT " + str(pf_cg.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("GONJUGATE_GRADIENT_" + str(pf_cg.count))
    #pf_save.rss_plot("NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)
    #pf.save_rss_plot_data(t_start, "NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)







 