
class pf_gd:

  h = 1.0e-5

  def run():
    
    # Start
    p = g.top_parameters[0][1]
    rss = pf_gd.rss(p)
    dp = pf_gd.df(p)
    dp = dp / max(dp)

    gamma = 1.0
    pf_gd.line_search(p, dp, gamma, rss)


    exit()

  def rss(p, top_parameters=False):
    pf_potential.update(p)
    return pf.get_rss(top_parameters)

  # Central Difference
  def df(p):
    dp = numpy.zeros((len(p),),)
    for i in range(len(p)):
      p_f = numpy.copy(p)
      p_b = numpy.copy(p)

      p_f[i] = p_f[i] + pf_gd.h
      p_b[i] = p_b[i] - pf_gd.h

      rss_f = pf_gd.rss(p_f)
      rss_b = pf_gd.rss(p_b)
      print(i, rss_b, rss_f)

      dp[i] = (pf_gd.rss(p_f) - pf_gd.rss(p_b)) / (2 * pf_gd.h)
    return dp

  
  def line_search(p, dp, gamma, best_rss):    
    gamma = 1.0
    while(gamma > 1.0e-24):
      p_test = numpy.copy(p)
      p_test = p_test - gamma * dp
      rss = pf_gd.rss(p_test)
      gamma = 0.25 * gamma
      print(gamma, rss, best_rss)
      if(rss < best_rss):
        p_best = p_test
        best_rss = rss
        best_gamma = gamma
        

