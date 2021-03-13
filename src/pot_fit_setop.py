class pf_setop:

  def run():

    for n in range(g.fit['start_random']):
      p = pf_parameters.random_p()
      pf_potential.update(p)
      rss = pf.get_rss()

    step_enhance_steps = 10
    
    # Loop through entire top list
    p_enhance = []
    for i in range(len(g.top_100)):
      p_enhance.append(g.top_100[i][1])

    # Make new list of top parameters
    new_list = pf_setop.process_list(p_enhance, step_enhance_steps)

    # Empty and re-populate top
    g.top_100 = []
    for i in range(len(new_list)):
      pf_potential.update(new_list[i])
      rss = pf.get_rss()



       

  
  def process_list(p_list, steps = 100):

    p_new_list = []

    for n in range(len(p_list)):

      g.pfdata['stage'] = 'STEP_' + str(n+1)

      p_current = numpy.copy(p_list[n])
      rss_current = pf.get_rss()
      rss_start = rss_current

      for n in range(steps):
        p_new = numpy.copy(p_current)
        p_new = p_new + pf_setop.random_step(0.0001)

        pf_potential.update(p_new)
        rss = pf.get_rss(False)

        if(rss < rss_current):
          p_current = p_new
          rss_current = rss
      
      p_new_list.append(p_current)

    return p_new_list


  def random_step(m = 0.01):
    # Get random parameters 0 to 1
    return m * (0.5 - numpy.random.rand(g.pfdata['params']['count']))