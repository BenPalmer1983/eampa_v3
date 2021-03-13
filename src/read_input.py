
from globals import g





class read_input:
 
 
 
  def run():
    main.log_title("Read Input")
  
    read_input.run_type()
    read_input.wd()
    read_input.rss_weights()
    read_input.rss_max_density()
    read_input.fitting()
    read_input.fit()
    read_input.fit_results()
    read_input.bp()
    read_input.mask()
    read_input.dft_ea()


  # READ TYPE
  def run_type():
  
    # DEFAULT
    g.run_type = 'efs'
    
    # TRY READING
    try:
      g.run_type = g.inp['run']['type'].lower().strip()
    except:
      pass
      
    # SAVE
    main.log(g.run_type)



  # READ TYPE
  def wd():
    g.wd_type = {
    'option': 1,
    }    
    # TRY READING
    for k in g.wd_type.keys():
      try:
        g.wd_type[k] = float(g.inp['wd'][k])
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
    'negec': 1.0,
    }
    
    # TRY READING
    for k in g.rss_weights.keys():
      try:
        g.rss_weights[k] = float(g.inp['rss_weights'][k])
      except:
        pass
      # READ 
      
      
    # SAVE
    main.log(std.dict_to_str(g.rss_weights))
      
      
      
   
  # READ 
  def rss_max_density():   
      
    # DEFAULT
    g.rss_max_density = {          
    'min': 0.2,
    'max': 0.8,
    'scale_factor': 10.0,
    'scale_exponent': 4.0,
    'zero_density_factor': 1.0e8,
    }
    
    # TRY READING
    for k in g.rss_max_density.keys():
      try:
        g.rss_max_density[k] = float(g.inp['rss_max_density'][k])
      except:
        pass
      # READ 
      
      
    # SAVE
    main.log(std.dict_to_str(g.rss_max_density))
    
      
   
  # READ 
  def fitting():   
      
    # DEFAULT
    g.fitting = {
    'oversized_parameters': [10.0,0.05,0.05,0.05],   
    'top_parameters': 100,
    'load_top_parameters': 10,
    }
    if('fitting' in g.inp.keys()):
      for k in g.fitting.keys():
        try:
          g.fitting[k] = g.inp['fitting'][k]
        except:
          pass

      
  def fit():
  
    # List of fit types
    g.fit = []
    if('fit' in g.inp.keys()):
      g.fit.append(read_input.read_fit('fit'))
    for i in range(100):
      if("fit" + str(i) in g.inp.keys()):
        g.fit.append(read_input.read_fit("fit" + str(i)))

  def read_fit(inp_key):  
    
    fit_data = {
    'type': None,
    'random_size': 100,
    'cycles': 0,
    'gens': 0,
    'spline_cycles': 0,
    'spline_gens': 0,
    'pop_size': 20,
    'fresh_size': 10,
    'gen_variation_multiplier': 0.1,
    'exct_factor': 0.5,
    'exct_every': 5,
    'exct_var': 0.1,
    'exct_top_bias': 0.5,
    'rescale_density': 0,
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
    'sa_loops_t': 10,
    'sa_loops_i': 100,
    'sa_temp_start': 10,
    'sa_temp_end': 0.1,
    'sa_step': 0.01,
    'sa_step_factor': 0.3,
    }

    # TRY READING
    for k in fit_data.keys():
      try:
        fit_data[k] = g.inp[inp_key][k]
      except:
        pass
        
    # POP SIZE - must be even
    if(fit_data['pop_size'] < 2):
      fit_data['pop_size'] = 2
    if(fit_data['pop_size'] % 2 != 0):
      fit_data['pop_size'] = fit_data['pop_size'] + 1
      
    # FRESH SIZE - must be even
    if(fit_data['fresh_size'] < 2):
      fit_data['fresh_size'] = 2
    if(fit_data['fresh_size'] % 2 != 0):
      fit_data['fresh_size'] = fit_data['fresh_size'] + 1

    if(fit_data['type'] is None):
      fit_data['type'] = 'RANDOM'


    return fit_data



      
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
      
    # SAVE
    main.log(std.dict_to_str(g.fit_results))
    
    

  def bp():
      
    # DEFAULT
    g.bp_input = {
    'dir': '',
    'bp_file': None,
    'eos_size': 10,
    'eos_strain': 0.005,
    'ec_size': 10,
    'ec_strain': 0.005,
    }

    # TRY READING
    for k in g.bp_input.keys():
      try:
        g.bp_input[k] = g.inp['bp'][k]
      except:
        pass
        
    # SAVE
    main.log(std.dict_to_str(g.bp_input))


  
  def mask():
    g.mask = {}
    if('mask' in g.inp):
      try:
        for k in g.inp['mask']:
          g.mask[k.upper()] = g.inp['mask'][k].upper()
      except:
        pass
        
    # SAVE
    main.log(std.dict_to_str(g.mask))




  def dft_ea():
  
    dftea = {}
    if('dft' in g.inp):
      try:
        for k in g.inp['dft']:
          dftea[k.upper()] = g.inp['dft'][k]
      except:
        pass
        
    
    for k in dftea.keys():
      label_str, label_id = labels.add(k)  
      
      atom_count = int(dftea[k][0])
      relaxed_energy = float(dftea[k][1])
      relaxed_energy_unit = str(dftea[k][2])
      coh_energy = float(dftea[k][3])
      coh_energy_unit = str(dftea[k][4])
      
      relaxed_dft_ev = units.convert(relaxed_energy_unit, "EV", relaxed_energy / atom_count) # relaxed per atom in eV
      coh_ev = units.convert(coh_energy_unit, "EV", coh_energy)
      apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
      g.dft_energy_adjustments[label_id] = {
                                            'label_id': label_id,
                                            'label_text': label_str,
                                            'atom_count': atom_count,
                                            'relaxed_energy': relaxed_energy,
                                            'relaxed_energy_unit': relaxed_energy_unit,
                                            'coh_energy': coh_energy,
                                            'coh_energy_unit': coh_energy_unit,
                                            'calc_relaxed_dft_ev': relaxed_dft_ev,
                                            'calc_coh_ev': coh_ev,
                                            'calc_apaev': apaev,
                                           }
      
    # Save  
    main.log(std.dict_to_str(g.dft_energy_adjustments))
  
  
  
  
  
  
  
  