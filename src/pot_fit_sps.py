"""
Simple single parameter step
No direction
Temperature loop from start to end temperature
Step size decreases with temperature

"""

class pf_sps:

  t = None
  count = 0

  def run(fit):
    for repeat in range(fit['repeat']):
      pf_sps.run_sps(fit)

  def run_sps(fit):
    pf_sps.count = pf_sps.count + 1

    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 

    g.pfdata['stage'] = 'Single Parameter Step ' + str(pf_sps.count)
    g.pfdata['stage_brief'] = 'SPS' + str(pf_sps.count)

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

    f = 0.01
    for n_out in range(5):
      for n_in in range(10):
        p_var = potential.get_p_var()
        for pn in range(len(p_var)):
          if(p_var[pn] != 0.0):
            p_new = numpy.copy(p_best)
            p_new[pn] = p_new[pn] + f * p_var[pn] * (0.5 - numpy.random.rand())
            potential.update(p_new)
            rss_new = pf.get_rss(False)
            if(rss_new < rss_best):
              rss_best = rss_new
              p_best = numpy.copy(p_new)
              rss_plot.append([time.time()-t_start, rss_best])
      f = f * 0.1
    rss_plot.append([time.time()-t_start, rss_best])

    pf.summary_line("SINGLE_PARAMETER_STEP " + str(pf_sps.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("SINGLE_PARAMETER_STEP" + str(pf_sps.count))
    pf_save.rss_plot("SINGLE_PARAMETER_STEP_" + str(pf_sps.count), rss_plot)
    pf.save_rss_plot_data(t_start, "SINGLE_PARAMETER_STEP_" + str(pf_sps.count), rss_plot)

    exit()




