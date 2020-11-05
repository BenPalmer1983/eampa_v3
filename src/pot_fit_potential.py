


class pf_potential:

  def update(p, no_rescale=False):    
    # Store
    g.pfdata['params']['current'][:] = p[:]
  
    # Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
        # Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size']  
        
        g.pot_functions['functions'][fn]['s_nodes'][:,1] = p[a:b]
        potential.make_spline_points_inner(fn)
        
        # Make Analytic Points
        #g.pot_functions['functions'][fn]['s_params'][:] = p[a:b]
        # LOAD ORIGINAL
        #g.pot_functions['functions'][fn]['points'] = numpy.copy(g.pot_functions['functions'][fn]['points_original'])   
        # VARY SPLINE
        #potential.vary_tabulated_points(fn, p[a:b])
        
        # Update a
        a = b   
        
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        # Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        # Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)
        a = b    
    

    # Rescale density functions 
    if(g.fit['rescale_density'] == 2 and no_rescale == False):
      rescale_density.run()    
    
    # Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    


  def take_density(p, p_dens): 
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      b = a + g.pot_functions['functions'][fn]['fit_size'] 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE       
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          p[a:b] = copy.deepcopy(p_dens[a:b])
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC          
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          p[a:b] = copy.deepcopy(p_dens[a:b])
      a = b
    return p