class pf_enhance:


  def run():
  
  
    if(g.pfdata['enhance']['top'] == 0):
      return 0
  
    g.pfdata['enhance']['counter'] += 1
    if(g.pfdata['enhance']['every'] == g.pfdata['enhance']['counter']):    
      g.pfdata['enhance']['event_count'] += 1
      g.pfdata['enhance']['counter'] = 0
    else:
      return 0
  
    pop_count = g.pfdata['enhance']['top']
    if(pop_count > g.fit['pop_size']):
      pop_count = g.fit['pop_size']
      
      
    for p in range(pop_count):
      print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
    
    for p in range(pop_count):
      print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
      params = gd.opt(pf_enhance.gd_rss, g.pfdata['params']['pop'][p,:])
      if(gd.rss_out < g.pfdata['params']['pop_rss'][p]):
        g.pfdata['params']['pop_rss'][p] = gd.rss_out
        g.pfdata['params']['pop'][p,:] = numpy.copy(params)
  
  
    #enhance_freq     #enhance_freq 
    
  def gd_rss(params):
    pf_potential.update(params)
    rss = pf.get_rss()
    return rss 