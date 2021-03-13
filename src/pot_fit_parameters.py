
class pf_parameters:

            

  def random_p(c=0.0, m=None):
  
    # Multiply the range
    if(m == None):
      m = 1.0

    # Randomly pick over sized parameters (defined in input file)
    rn = numpy.random.uniform()    
    b = g.fitting['oversized_parameters'][0]
    prb = 1.0
    for i in range(1, len(g.fitting['oversized_parameters']),1):
      prb = prb * g.fitting['oversized_parameters'][i]
      if(rn<prb):
        m = m * b
      else:
        break

    # Get the upper and lower range
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower
    
    # If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      center = lower + 0.5 * p_range 
    else:
      center = numpy.copy(c)
    
    # Multiply range
    m_range = m * p_range
    
    # Get random parameters 0 to 1
    r = numpy.random.rand(p_count)
    
    # New parameters
    p_new = center + (r - 0.5) * m_range
    
    # Return
    return p_new
    
    
    
  def random_variation(p_in, variation):
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower

    # p_out
    p_out = numpy.copy(p_in)
    p_out = p_out +  variation * p_range * (numpy.random.rand(p_count) - 0.5)

    return p_out


  def mutate(p_in, variation):
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower

    # p_out
    p_out = numpy.copy(p_in)

    












###############################################################################################
#
# Delete some day
#
###############################################################################################    

    
    #e = 0
    #if(g.pfdata['generation']['counter']>0):
    #  e = g.pfdata['generation']['counter'] - 1
    #p_new = c + (r - 0.5) * m_range * (g.fit['gen_var_factor'])**(e)
    
    
    

  def get_p_from_pool():
    pmax = g.pfdata['pool']['pmax']
    pn = g.pfdata['pool']['pn'] % pmax
    p = copy.deepcopy(g.pfdata['pool']['params'][pn, :])                     # Copy parameters
    g.pfdata['pool']['params'][pn, :] = pf_parameters.random_p(p, 0.01)      # Disturb
    g.pfdata['pool']['pn'] = g.pfdata['pool']['pn'] + 1                      # Increment
    return p                                                                 # Return
  
  
  
  def setup_pool():
    width = potential.parameter_count()
    g.pfdata['pool']['sane_a'] = g.fit['sane_seeds_a']
    g.pfdata['pool']['sane_b'] = g.fit['sane_seeds_b']
    g.pfdata['pool']['pmax'] = g.fit['pool_size']
    
    if(g.pfdata['pool']['pmax'] < (g.pfdata['pool']['sane_a'] + g.pfdata['pool']['sane_b'])):
      g.pfdata['pool']['pmax'] = g.pfdata['pool']['sane_a'] + g.pfdata['pool']['sane_b']
    g.pfdata['pool']['params'] = numpy.zeros((g.pfdata['pool']['pmax'],width,),)
  
  
  def make_pool(p = None):
    print("Make Pool")
    main.log_title("Make Pool")    
    main.log("pool size:     " + str(g.fit['pool_size']))
    main.log("pop size:      " + str(g.fit['pop_size']))
    main.log("fresh size:    " + str(g.fit['fresh_size']))
    
    
    if(p == None):
      if(g.pfdata['generation']['counter'] == 0):
        c = 0.0
        main.log("c:    " + str(c))
      else:
        c = g.pfdata['params']['best']
        
        str_c = ''
        for i in range(len(c)):
          str_c = str_c + str(c[i]) + ' '
        main.log("c:    " + str(str_c))
        
   
   
  
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    pmax = g.fit['pool_size']
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    sa = g.pfdata['pool']['sane_a']
    sb = g.pfdata['pool']['sane_b']
    width = potential.parameter_count()
    #width = potential.parameter_count()  
    #g.pfdata['pool']['params'] = numpy.zeros((pmax,width,),)
    #g.pfdata['pool']['pmax'] = pmax
  
  
    # Any Density
    if(g.fit['rescale_density'] == 0):
      print("Creating Pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(c, w)
        g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        #print(pn)
    elif(g.fit['rescale_density'] == 1 or g.fit['rescale_density'] == 2):
    
      sane_a = numpy.zeros((sa,width,),)
      sane_b = numpy.zeros((sb,width,),)
      
      print("Make Sane A")
      main.log("make sane set A")
      pn = 0
      while(pn < sa):
        sane_a[pn,:] = pf_cycle.random_p(c)
        a = 0
        for fn in range(len(g.pot_functions['functions'])): 
          b = a + g.pot_functions['functions'][fn]['fit_size'] 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            loop = True
            while(loop):
              pf_potential.update(sane_a[pn,:], True)
              rho = pf_parameters.estimate_density(fn)
              if( rho > 0.0 and rho <=1.0):
                loop = False
              else:
                ptemp = pf_cycle.random_p(c) 
                sane_a[pn,a:b] = numpy.copy(ptemp[a:b])
          a = b
        pn = pn + 1
              
      main.log(str(pn))      
      print("Make Sane B")
      main.log("make sane set B")
      pn = 0        
      while(pn < sb):
        loop = True
        while(loop):
          loop = False
          r = numpy.random.rand(width)
          params = (1.0 + 0.1 * (0.5-r)) * sane_a[pn%sa,:]
          pf_potential.update(params, True)
          for fn in range(len(g.pot_functions['functions'])): 
            if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
              rho = pf_parameters.estimate_density(fn)
              if(not(rho > 0.0 and rho <=1.0)): 
                loop = True             
          sane_b[pn,:] = params
        pn = pn + 1
      main.log(str(pn))
        
      
      print("Make Pool")
      main.log("make pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - (1+sa+sb))
      pn = 0
      while(pn<pmax):
        p = pf_parameters.random_p(c, w)
        if(g.fit['sane_fraction']>numpy.random.rand()):
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(pf_potential.take_density(p, sane_b[pn%sb,:]))  
        else:
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        
      #numpy.savetxt('test.out', g.pfdata['pool']['params'], delimiter=',') 
      #exit()  
        
    # Shuffle    
    print("Shuffle Pool")
    main.log("shuffle pool")
    numpy.random.shuffle(g.pfdata['pool']['params'])

    
    
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