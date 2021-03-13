


class pf_potential:

  def update(p, no_rescale=False):    
  
    # Update potential
    a = 0

    for fn in range(len(g.pot_functions['functions'])): 
      # Calc b
      b = a + g.pot_functions['functions'][fn]['fit_size']  
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE           
        g.pot_functions['functions'][fn]['s_nodes'][:,1] = p[a:b]
        potential.make_spline_points_inner(fn)       
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        # Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)

      # Update a
      a = b     
        

    # Rescale embedding data points to cover density range
    potential.rescale_embedding()

    # Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    
    # Save parameters
    potential.parameters = numpy.copy(potential.get_parameters())
    g.pfdata['p']['current'] = numpy.copy(potential.get_parameters())
 


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











"""
  
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
"""