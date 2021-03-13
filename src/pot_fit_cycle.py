
class pf_cycle:



    
  
  def run():  
  
    # Increment Cycle Counter
    g.pfdata['cycle']['counter'] += 1
    
    # Load from globals
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    
    #############################
    # Initialise first half
    # Search within ranges
    #############################
    
    g.pfdata['stage'] = 'Initialising Population'
    
    
    # Add top parameters (if there are any)
    
    
    cycle_start_parameters = []
    new_p = copy.deepcopy(g.pfdata['params']['start'][:])
    cycle_start_parameters.append(new_p)
    
    for i in range(len(g.top_parameters)):
      new_p = copy.deepcopy(g.top_parameters[i][1][:])
      cycle_start_parameters.append(new_p)

    
    
    n = 0
    for pn in range(pop_size):
      loop = True
      while(loop): 
        n = n + 1
        if(n <= len(cycle_start_parameters)):
          g.pfdata['params']['pop'][pn, :] = cycle_start_parameters[n-1][:]
        else:
          g.pfdata['params']['pop'][pn, :] = pf_parameters.random_p()
        
        # Try - if it fails or rss == None, try next
        pf_potential.update(g.pfdata['params']['pop'][pn, :])
        rss = pf.get_rss()
        try:
          # Update
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['pop_rss'][pn] = rss
        except:
          pass

      
 
            
    ###################################
    # LOOP THROUGH GENERATIONS
    ###################################
  
    for gen in range(g.fit['gens']):
      g.pfdata['generation']['counter'] += 1
      g.pfdata['stage'] = 'Loop Through Generations - Gen ' + str(g.pfdata['generation']['counter']) 
      
      # Run through a generation
      pf_generation.run()
      
      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  def random_p(c=0.0, m=1.0):
    # 
    p_count = g.pfdata['params']['count']
    lower = g.pfdata['params']['var'][0,:]
    upper = g.pfdata['params']['var'][1,:]
    range = upper - lower
    
    # If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      c = lower + 0.5 * range 
    
    # Multiply range
    m_range = m * range
    
    # Get random parameters 0 to 1
    r = numpy.random.rand(p_count)
    
    # New parameters
    e = 0
    if(g.pfdata['generation']['counter']>0):
      e = g.pfdata['generation']['counter'] - 1
    
    p_new = c + (r - 0.5) * m_range * (g.fit['gen_var_factor'])**(e)
    
    # Return
    return p_new
    
    
  def run_old():  
  
    # Increment Cycle Counter
    g.pfdata['cycle']['counter'] += 1
    
    # Load from globals
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    
    #############################
    # Initialise first half
    # Search within ranges
    #############################
    
    g.pfdata['stage'] = 'Initialising Population - First Half'
    for p in range(pop_size // 2):
      if(p == 0):
        g.pfdata['params']['pop'][p, :] = g.pfdata['params']['start'][:]
        pf_potential.update(g.pfdata['params']['pop'][p, :])
        rss = pf.get_rss()
        g.pfdata['params']['pop_rss'][p] = rss
      else:
        loop = True
        while(loop): 
          g.pfdata['params']['pop'][p, :] = pf_cycle.random_p(0.0, 1.0)
          try:
            pf_potential.update(g.pfdata['params']['pop'][p, :])
            rss = pf.get_rss()
            if(g.pfdata['rss']['current'] is not None):
              loop = False
              g.pfdata['params']['pop_rss'][p] = rss
          except:
            pass

        
    #############################
    # Initialise second half
    # Search a wider area
    #############################
  
    g.pfdata['stage'] = 'Initialising Population - Second Half'
    w = w_start
    w_inc = (w_end - w_start) / (pop_size // 2 - 1)
    for p in range(pop_size // 2, pop_size):
      w = w + w_inc
      loop = True
      while(loop): 
        g.pfdata['params']['pop'][p, :] = pf_cycle.random_p(0.0, w)
        try:
          pf_potential.update(g.pfdata['params']['pop'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['pop_rss'][p] = rss
        except:
          pass
            
    ###################################
    # LOOP THROUGH GENERATIONS
    ###################################
  
    for gen in range(g.fit['gens']):
      g.pfdata['generation']['counter'] += 1
      g.pfdata['stage'] = 'Loop Through Generations - Gen ' + str(g.pfdata['generation']['counter']) 
      
      # Run through a generation
      pf_generation.run()
      
      
      
      
      