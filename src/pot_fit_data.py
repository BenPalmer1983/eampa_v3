



class pf_data:
  
  def make():

    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    width = potential.parameter_count()
    
    top_size = 10
    if(top_size < 10):
      top_size = 10
    top = {'filled': False, 'counter': 0, 'size': top_size, 'rss': numpy.zeros((top_size,),), 'p': numpy.zeros((top_size, width,),),}
    top['rss'][:] = -1.0

    rss = { 'start': None,
            'current': None,
            'best': None,
            'counter': 0,
            'counter_successful': 0,
            'since_improvement': 0,
    }
    d = {
        'stage': 'Not Set', 
        'params': {'count': width, 
                   'start': numpy.zeros((width,),),
                   'var': numpy.zeros((2,width,),),
                   'best': numpy.zeros((width,),),
                   'pop': numpy.zeros((pop_size, width,),),
                   'pop_rss': numpy.zeros((pop_size,),),
                   'fresh': numpy.zeros((fresh_size, width,),),
                   'fresh_rss': numpy.zeros((fresh_size,),),
                   'children': numpy.zeros((pop_size + 2 * fresh_size, width,),),
                   'children_rss': numpy.zeros((pop_size + 2 * fresh_size,),),
                   'current': numpy.zeros((width,),),
                   'pop_merged': numpy.zeros((2 * pop_size + 3 * fresh_size, width,),),
                   'pop_merged_rss': numpy.zeros((2 * pop_size + 3 * fresh_size,),),
                  },
        'rss': rss,
        'cycle': {'counter': 0, 'total_cycles': g.fit['cycles'], 'spline_counter': 0,},
        'generation': {'counter': 0, 'total_generations': g.fit['gens'], 'spline_counter': 0,},
        'top': top,
        'extinction': {'event_count': 0, 'every': g.fit['exct_every'], 'counter': 0, 'top_bias': g.fit['exct_top_bias'],},
        'enhance': {'event_count': 0, 'every': g.fit['enhance_every'], 'top': g.fit['enhance_top'], 'counter': 0,},
        'best_hash': None,
        'bp': {},
        'efs': {},
        'bp_known': g.bp_known,
        'efs_known': g.efs_known,
        'bp_best': None,
        'time': {'start': g.times['start'], 'end': None},
    }
    
    
    # SAVE Starting Parameters
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        d['params']['start'][a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        d['params']['start'][a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b    
    
    
    # SAVE Variation 
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        d['params']['var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        d['params']['var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        d['params']['var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        d['params']['var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  
    
    return d





