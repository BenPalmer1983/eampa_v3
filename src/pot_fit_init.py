class pf_init:

  def run():  

    # If a spline fit, convert into a spline with the set number of nodes
    for fn in range(len(g.pot_functions['functions'])):       
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     
        potential.vary_tabulated_points(fn)

  
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    efs_calc.set_weights()             # Set weightings
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
    # Setup BP    
    bp_calc.init()
    bp_calc.set_weights()
    potential.bp_add_potentials()
    b_props.bp_add()
    bp_calc.get_known()

    g.pfdata = {}
    
    g.pfdata['stage'] = ''
    g.pfdata['last_stage'] = None

    # Set up RSS
    g.pfdata['rss'] = {'start': None,
                       'current': None,
                       'best': None,
                       'counter': 0,
                       'counter_successful': 0,
                       'since_improvement': 0,}


    g.pfdata['psize'] = potential.parameter_count()
    g.pfdata['p'] = {'current': None, 'best': None,}

    g.pfdata['bp'] = {'current': None, 'best': None,}

    g.pfdata['max_density'] = {'bp_current': None, 'efs_current': None, 'bp_best': None, 'efs_best': None,}
    


    # SAVE Starting Parameters
    g.pfdata['params_start'] = numpy.zeros((g.pfdata['psize'],),)
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        g.pfdata['params_start'][a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        g.pfdata['params_start'][a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b    
    
    
    # SAVE Variation 
    g.pfdata['params_var'] = numpy.zeros((2,g.pfdata['psize'],),)
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        g.pfdata['params_var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        g.pfdata['params_var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        g.pfdata['params_var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        g.pfdata['params_var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  


















