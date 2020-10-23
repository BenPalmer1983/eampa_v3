
from globals import g





class read_input:
 
 
 
  def run():
    read_input.run_type()
    read_input.rss_weights()
    read_input.fit()
    read_input.fit_results()



  # READ TYPE
  def run_type():
  
    # DEFAULT
    g.run_type = 'efs'
    
    # TRY READING
    try:
      g.run_type = g.inp['run']['type'].lower().strip()
    except:
      pass



  # READ 
  def rss_weights():
  
    # DEFAULT
    g.rss_weights = {
    'config': 1.0,
    'force': 1.0,
    'energy': 1.0,
    'stress': 1.0,
    'bp': 1.0,
    'a0': 1.0,
    'e0': 1.0,
    'b0': 1.0,
    'ec': 1.0,
    'g': 1.0,
    'e': 1.0,
    'v': 1.0,
    }
    
    # TRY READING
    for k in g.rss_weights.keys():
      try:
        g.rss_weights[k] = float(g.inp['rss_weights'][k])
      except:
        pass
      # READ 
      
      
      
      
  def fit():
  
    # DEFAULT
    g.fit = {
    'cycles': 1,
    'gens': 4,
    'spline_cycles': 0,
    'spline_gens': 0,
    'pop_size': 20,
    'fresh_size': 10,
    'exct_factor': 0.5,
    'exct_every': 5,
    'exct_var': 0.1,
    'exct_top_bias': 0.5,
    'rescale_density': 0,
    'rescale_min': 0.3,
    'rescale_max': 0.9,
    'rescale_default': 0.6,
    'wide_start': 0.5,
    'wide_end': 10.0,
    'mutate_chance': 0.01,
    'mutate_scale': 1.0,
    'no_clones': True,
    'no_clone_var': 0.05,
    'fresh_ws': 0.1,
    'fresh_we': 1.5,
    'enhance_every': 10,
    'enhance_top': 5,
    'gen_var_factor': 1.0,
    'pool_size': 1000,
    'sane_seeds_a': 50,
    'sane_seeds_b': 200,
    'sane_fraction': 0.5,
    }

    # TRY READING
    for k in g.fit.keys():
      try:
        g.fit[k] = g.inp['fit'][k]
      except:
        pass
        
    # POP SIZE - must be even
    if(g.fit['pop_size'] < 2):
      g.fit['pop_size'] = 2
    if(g.fit['pop_size'] % 2 != 0):
      g.fit['pop_size'] = g.fit['pop_size'] + 1
      
    # FRESH SIZE - must be even
    if(g.fit['fresh_size'] < 2):
      g.fit['fresh_size'] = 2
    if(g.fit['fresh_size'] % 2 != 0):
      g.fit['fresh_size'] = g.fit['fresh_size'] + 1
      
      
      
  def fit_results():
      
    # DEFAULT
    g.fit_results = {
    'results_dir': 'results',
    'pot_dir': None,
    'pot_name': None,
    }
    
    # TRY READING
    for k in g.fit_results.keys():
      try:
        g.fit_results[k] = g.inp['fit_results'][k]
      except:
        pass
      
