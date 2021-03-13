

class pf_steps:


  
  def run():

    p = numpy.copy(g.pfdata['params']['start'])
    rss = pf.get_rss()
    
    for n in range(5):
      p_new = numpy.copy(p)
      p_new = p_new + pf_steps.random_step()

      pf_potential.update(p_new)
      rss_new = pf.get_rss()
  
      if(rss_new < rss):
        p = numpy.copy(p_new)
        rss_new = rss


  def random_step(m = 0.01):
    # Get random parameters 0 to 1
    return m * (0.5 - numpy.random.rand(g.pfdata['params']['count']))
    
    

  
  def process_list(p_list, steps = 100):

    p_new_list = []

    for n in range(len(p_list)):

      g.pfdata['stage'] = 'STEP_' + str(n+1)

      p_current = numpy.copy(p_list[n])
      rss_current = pf.get_rss()
      rss_start = rss_current

      for n in range(steps):
        p_new = numpy.copy(p_current)
        p_new = p_new + pf_steps.random_step(0.0001)

        pf_potential.update(p_new)
        rss = pf.get_rss(False)

        if(rss < rss_current):
          p_current = p_new
          rss_current = rss
      
      p_new_list.append(p_current)

    return p_new_list
