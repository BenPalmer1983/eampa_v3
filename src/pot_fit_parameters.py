
class pf_parameters:


  def get_p():
    pmax = g.pfdata['pool']['pmax']
    pn = g.pfdata['pool']['pn'] % pmax
    p = copy.deepcopy(g.pfdata['pool']['params'][pn, :])                     # Copy parameters
    g.pfdata['pool']['params'][pn, :] = pf_parameters.random_p(p, 0.01)      # Disturb
    g.pfdata['pool']['pn'] = g.pfdata['pool']['pn'] + 1                      # Increment
    return p                                                                 # Return
  
  def make_pool():
    print("Make Pool")
  
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    pmax = g.fit['pool_size']
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    width = potential.parameter_count()
  
    g.pfdata['pool']['params'] = numpy.zeros((pmax,width,),)
    g.pfdata['pool']['pmax'] = pmax
  
  
    # Any Density
    if(g.fit['rescale_density'] == 0):
      print("Creating Pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(0.0, w)
        g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        #print(pn)
    elif(g.fit['rescale_density'] == 1 or g.fit['rescale_density'] == 2):
      sa = g.fit['sane_seeds_a']
      sb = g.fit['sane_seeds_b']
    
      sane_a = numpy.zeros((sa,width,),)
      sane_b = numpy.zeros((sb,width,),)
      
      print("Make Sane A")
      pn = 0
      while(pn < sa):
        params = pf_cycle.random_p()
        pf_potential.update(params, True)
        for fn in range(len(g.pot_functions['functions'])): 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            rho = pf_parameters.estimate_density(fn)
            if( rho > 0.0 and rho <=1.0):
              sane_a[pn,:] = params
              pn = pn + 1
              
      print("Make Sane B")
      pn = 0        
      while(pn < sb):
        r = numpy.random.rand(width)
        params = (1.0 + 0.1 * (0.5-r)) * sane_a[pn%sa,:]
        pf_potential.update(params, True)
        for fn in range(len(g.pot_functions['functions'])): 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            rho = pf_parameters.estimate_density(fn)
            if( rho > 0.0 and rho <=1.0):
              sane_b[pn,:] = params
              pn = pn + 1
        
      
      print("Make Pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(0.0, w)
        if(g.fit['sane_fraction']>numpy.random.rand()):
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(pf_potential.take_density(p, sane_b[pn%sb,:]))  
        else:
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        
    # Shuffle    
    print("Shuffle Pool")
    numpy.random.shuffle(g.pfdata['pool']['params'])

            

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
    
    
  def estimate_density(fn):  
    #print(g.pot_functions['functions'][fn]['points'][:,1])
    r = numpy.zeros((7,),)
    rn = numpy.zeros((7,),)
    r[0] = 7.48332e0
    r[1] = 6.32456e0
    r[2] = 6.92821e0
    r[3] = 4.89898e0
    r[4] = 5.65686e0
    r[5] = 2.82843e0
    r[6] = 4.0e0
    rn[0] = 48
    rn[1] = 24
    rn[2] = 8
    rn[3] = 24
    rn[4] = 12
    rn[5] = 12
    rn[6] = 6    
    rho = 0.0    
    for i in range(7):
      y = interp.search_x(r[i], g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
      rho = rho + rn[i] * y
      #print(rho)
    return rho
    
    
  def p_diff():
    return (g.pfdata['params']['var'][1,:] - g.pfdata['params']['var'][0,:])