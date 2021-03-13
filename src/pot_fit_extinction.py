

class pf_extinction:




  def run():
  
    g.pfdata['extinction']['counter'] += 1
    if(g.pfdata['extinction']['every'] == g.pfdata['extinction']['counter']):    
      g.pfdata['extinction']['event_count'] += 1
      g.pfdata['extinction']['counter'] = 0
    else:
      return 0
  
    pf_extinction.best = g.pfdata['params']['pop_merged_rss'][0]
    pf_extinction.worst = g.pfdata['params']['pop_merged_rss'][-1]
    
    top_ten = copy.deepcopy(g.pfdata['params']['pop_merged'][0:10])
        
    #print("Extinction")
    #print(pf_extinction.best, pf_extinction.worst)
    
    # Cause extinction of poorly fitting parameters
    for i in range(len(g.pfdata['params']['pop_merged_rss'][:])):
      chance = pf_extinction.chance(g.pfdata['params']['pop_merged_rss'][i])
      r = numpy.random.uniform(0.0,1.0)
      if(r < chance):       
        #print("Change", i, g.pfdata['params']['pop_merged_rss'][i])
        loop = True
        while(loop):
          if(numpy.random.uniform(0.0,1.0)<=g.fit['exct_top_bias']):
            rn = 0
          else:
            rn = numpy.random.randint(0,9)
          new_p = pf_cycle.random_p(top_ten[rn], g.fit['exct_var'])
          pf_potential.update(new_p)
          rss = pf.get_rss()
          if(rss != None):
            #print(rss, g.pfdata['params']['pop_merged_rss'][i])
            if(rss < g.pfdata['params']['pop_merged_rss'][i]):
              g.pfdata['params']['pop_merged'][i,:] = new_p
              g.pfdata['params']['pop_merged_rss'][i] = rss
            loop = False

    # Sort
    pf_extinction.sort_merged()
    #print(g.pfdata['params']['pop_merged_rss'])

  def chance(rss):
    return g.fit['exct_factor'] * ((rss - pf_extinction.best) / (pf_extinction.worst - pf_extinction.best)) 


  def sort_merged():
    sort.sort_1d_dp_asc(g.pfdata['params']['pop_merged_rss'][:])
    g.pfdata['params']['pop_merged_rss'] = sort.apply_keytable_1d_dp(g.pfdata['params']['pop_merged_rss'])
    g.pfdata['params']['pop_merged'] = sort.apply_keytable_2d_dp(g.pfdata['params']['pop_merged'])


  