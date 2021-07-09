
# Basin Hopping

class pf_bh:

  count = 0
  p_full = None

  def run(fit):
    for repeat in range(fit['repeat']):
      pf_bh.run_bh(fit)

  def get_rss(ps):
    p = pf_bh.grow_p(ps)
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
        p_out.append(pf_bh.p_full[i])
    p_out = numpy.asarray(p_out)
    return p_out
    

  def run_bh(fit):
    pf_bh.count = pf_bh.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 


    g.pfdata['stage'] = 'Basin-Hopping ' + str(pf_bh.count)
    g.pfdata['stage_brief'] = 'BH' + str(pf_bh.count)
  
    pf_bh.p_full = numpy.copy(g.pfdata['p']['best'])
    p = numpy.copy(g.pfdata['p']['best'])

    # Shrink to non-fixed parameters
    ps = pf_bh.shrink_p(p)
    res = basinhopping(pf_bh.get_rss, ps, niter=float(fit['bh_niter']), T=1.0, stepsize=0.5)
    ps = numpy.copy(res['x'])
    p = pf_bh.grow_p(ps)
    potential.update(p)
    rss = pf.get_rss(False)


    pf.summary_line("BASIN HOPPING " + str(pf_bh.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("BASIN_HOPPING_" + str(pf_bh.count))
    #pf_save.rss_plot("NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)
    #pf.save_rss_plot_data(t_start, "NELDER_MEAD_" + str(pf_neldermead.count), rss_plot)







 