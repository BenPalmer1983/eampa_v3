class pf_densities:
  
  def run():
    sa = len(g.pfdata['sane']['seeds_a'][:,0])
    sb = len(g.pfdata['sane']['seeds_b'][:,0])
    
    
    p_initial = copy.deepcopy(g.pfdata['params']['current'])
    
    print("Making Seed for Pool")
    pn = 0
    while(pn < sa):
      params = pf_cycle.random_p()
      pf_potential.update(params)
      for fn in range(len(g.pot_functions['functions'])): 
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          rho = pf_densities.estimate_density(fn)
          if( rho >= 0.0 and rho <=1.2):
            g.pfdata['sane']['seeds_a'][pn,:] = params
            pn = pn + 1
       
    print("Mutating Pool Seed for Pool")     
    for n in range(sb):
      p = random_p(g.pfdata['sane']['seeds_a'][n%sa,:], m=0.01)
    #print(g.pfdata['sane']['seeds'])     
 
  
  

  
  
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