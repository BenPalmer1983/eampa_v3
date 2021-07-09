


class pf_potential:

  # Analytic potentials only
  def update(p, no_rescale=False):    
    
    # Estimate densities to rescale embedding functions
    rho = pf_potential.estimate_densities()
  
    # Change upper bound for embedding functions
    """
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 3):
        g.pot_functions['functions'][fn]['a_u'] = 1.0
        a = g.pot_functions['functions'][fn]['a']
        f_group = g.pot_functions['functions'][fn]['f_group']
        try:
          rho_max = rho[a][f_group] * potential.embedding_range_factor
          if(rho_max > 1.0):
            g.pot_functions['functions'][fn]['a_u'] = rho_max 
            #print(fn, rho_max)
          else:
            g.pot_functions['functions'][fn]['a_u'] = 1.0
        except:
          pass
    """

    # Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      # Calc b
      b = a + g.pot_functions['functions'][fn]['fit_size']  
      # Make Analytic Points
      g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
      potential.make_analytic_points_inner(fn)
      # Update a
      a = b     
        

    # Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    
    # Save parameters
    potential.parameters = numpy.copy(potential.get_parameters())
    g.pfdata['p']['current'] = numpy.copy(potential.get_parameters())
 

  def estimate_densities():
    out = {}
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        rho_max = rescale_density.estimate_density(fn)
        a = g.pot_functions['functions'][fn]['a']
        f_group = g.pot_functions['functions'][fn]['f_group']
        if(a not in out.keys()):
          out[a] = {}
        out[a][f_group] = rho_max
    return out

  def estimate_density(fn):
    # Symetric simple cubic used as test case
    sc = {0.86602539999999995: 8, 0.80039053000000004: 24, 0.75: 24, 0.71807032999999998: 24, 0.70710678000000005: 12, 0.72886899000000005: 24, 0.67314560000000001: 48, 0.63737743999999996: 48, 0.625: 24, 0.61237244000000002: 24, 0.57282195999999996: 48, 0.55901699000000005: 24, 0.53033008999999998: 36, 0.51538819999999996: 48, 0.5: 6, 0.64951904999999999: 8, 0.58630196999999995: 24, 0.54486237000000004: 24, 0.46770717000000001: 48, 0.45069390999999998: 24, 0.4145781: 24, 0.39528470999999998: 24, 0.375: 30, 0.43301269999999997: 8, 0.35355339000000002: 12, 0.30618622000000001: 24, 0.27950849999999999: 24, 0.25: 6, 0.21650634999999999: 8, 0.17677670000000001: 12, 0.125: 7}
    a_prim = 0.125
    a = [2.0,3.0,4.0,5.0] 
    m_rho = 0.0
    for a0 in a:
      a0 = (1 / a_prim) * a0
      rho = 0.0    
      for k in sc.keys():
        y = interp.search_x(a0 * k, g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
        rho = rho + sc[k] * y
      m_rho = max(m_rho, rho)
    return m_rho




  def rescale_parameters(p):
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      # Calc b
      b = a + g.pot_functions['functions'][fn]['fit_size']  
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        #a

        """
        # Symetric simple cubic used as test case
        sc = {0.86602539999999995: 8, 0.80039053000000004: 24, 0.75: 24, 0.71807032999999998: 24, 0.70710678000000005: 12, 0.72886899000000005: 24, 0.67314560000000001: 48, 0.63737743999999996: 48, 0.625: 24, 0.61237244000000002: 24, 0.57282195999999996: 48, 0.55901699000000005: 24, 0.53033008999999998: 36, 0.51538819999999996: 48, 0.5: 6, 0.64951904999999999: 8, 0.58630196999999995: 24, 0.54486237000000004: 24, 0.46770717000000001: 48, 0.45069390999999998: 24, 0.4145781: 24, 0.39528470999999998: 24, 0.375: 30, 0.43301269999999997: 8, 0.35355339000000002: 12, 0.30618622000000001: 24, 0.27950849999999999: 24, 0.25: 6, 0.21650634999999999: 8, 0.17677670000000001: 12, 0.125: 7}
        a_prim = 0.125
        a0_list = [2.0,3.0,4.0,5.0] 
        m_rho = 0.0
        for a0 in a0_list:
          a0 = (1 / a_prim) * a0
          rho = 0.0    
          for k in sc.keys():
            y = g.pot_functions['functions'][fn]['f'](a0 * k,  
                                                      p[a:b],  
                                                      g.pot_functions['functions'][fn]['a_params_fixed'])
            rho = rho + sc[k] * y
          m_rho = max(m_rho, rho)
        """ 
            
      # Update a
      a = b     









  """
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