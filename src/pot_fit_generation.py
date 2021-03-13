"""
Shuffle parents
Breed within population
Breed with fresh pop



"""


class pf_generation:

  parents = None
  pop_size = None
  pop_size_half = None
  fresh_size = None

  def run():
    main.log_hr()
    main.log("Generation " + str(g.pfdata['generation']['counter']))
    main.log_hr()

  
    pf_generation.pop_size = g.fit['pop_size']
    pf_generation.pop_size_half = g.fit['pop_size'] // 2
    pf_generation.fresh_size = g.fit['fresh_size']

    # Make array of parents and shuffle
    pf_generation.parents = numpy.arange(pf_generation.pop_size)
    numpy.random.shuffle(pf_generation.parents)
    
    c = 0
    # Breed Population
    g.pfdata['stage'] = 'Breed Population'  
    for p in range(pf_generation.pop_size_half):  
      loop = True
      while(loop):
        loop = pf_generation.breed_event(p, c, 'p+p')      
      c = c + 2
       
    # Make Fresh  
    g.pfdata['stage'] = 'Initialising Fresh Population'   
    pf_generation.make_fresh()
       
    # Reshuffle
    numpy.random.shuffle(pf_generation.parents)
    
    # Breed Population-Fresh 
    g.pfdata['stage'] = 'Breed With Fresh Population'  
    for p in range(pf_generation.fresh_size):  
      loop = True
      while(loop):
        loop = pf_generation.breed_event(p, c, 'p+f')      
      c = c + 2
      
    # Merge populations and select top to form next generation   
    g.pfdata['stage'] = 'Merge Populations'   
    pf_generation.merge()
    
    
    # Extinction Event
    g.pfdata['stage'] = 'Extinction'   
    pf_extinction.run()
    
    # Load best back into population
    g.pfdata['params']['pop'][:, :] = g.pfdata['params']['pop_merged'][0:pf_generation.pop_size, :]
    g.pfdata['params']['pop_rss'][:] = g.pfdata['params']['pop_merged_rss'][0:pf_generation.pop_size]
    
    # Enhance
    g.pfdata['stage'] = 'Enhance'   
    pf_enhance.run()


    #g.pfdata['generation']['counter']
    #pf_setop
    
    
    
    
    cstr = str(g.pfdata['generation']['counter'])
    while(len(cstr)<5):
      cstr = '0' + cstr    
    gen_dir = g.dirs['wd'] + '/fitting/genetic/' +  cstr
    std.make_dir(gen_dir)
  
    
    potential_output.full(gen_dir)
    
    
    
    
    
  def breed_event(p, c, opt):
    pc = g.pfdata['params']['count']
    
    pa = pf_generation.parents[p]
    if(opt == 'p+p'):  
      pb = pf_generation.parents[p + pf_generation.pop_size_half]
      pb_array = 'pop'
    if(opt == 'p+f'): 
      pb = p
      pb_array = 'fresh'
    
    
    ca = c
    cb = c + 1
    
    # Breed 
    state = pf_generation.get_state()
    for i in range(pc):
      state = pf_generation.get_state(state)
      if(state):
        g.pfdata['params']['children'][ca, i] = g.pfdata['params']['pop'][pa, i]
        g.pfdata['params']['children'][cb, i] = g.pfdata['params'][pb_array][pb, i]
      else:      
        g.pfdata['params']['children'][cb, i] = g.pfdata['params']['pop'][pa, i]
        g.pfdata['params']['children'][ca, i] = g.pfdata['params'][pb_array][pb, i]
      
    # Mutate
    g.pfdata['params']['children'][ca, :] = pf_generation.mutate(g.pfdata['params']['children'][ca, :], g.fit['mutate_chance'])
    g.pfdata['params']['children'][cb, :] = pf_generation.mutate(g.pfdata['params']['children'][cb, :], g.fit['mutate_chance'])
    
    # No clones - Child A
    for pn in range(pf_generation.pop_size):
      if((g.pfdata['params']['children'][ca, :] == g.pfdata['params']['pop'][pn,:]).all()):
        r = numpy.random.rand(pc)   
        g.pfdata['params']['children'][ca, :] = g.pfdata['params']['children'][ca, :] * g.fit['no_clone_var'] * (0.5 - r)
        break
        
    # No clones - Child B
    for pn in range(pf_generation.pop_size):
      if((g.pfdata['params']['children'][cb, :] == g.pfdata['params']['pop'][pn,:]).all()):
        r = numpy.random.rand(pc)   
        g.pfdata['params']['children'][cb, :] = g.pfdata['params']['children'][cb, :] * g.fit['no_clone_var'] * (0.5 - r)
        break


    loop = False 
    
    pf_potential.update(g.pfdata['params']['children'][ca, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
    pf_potential.update(g.pfdata['params']['children'][cb, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
     
    if(loop == False):
      g.pfdata['params']['children_rss'][ca] = rss
      g.pfdata['params']['children_rss'][cb] = rss
      
      
    # Loop again or not
    return loop
    
    
    
  def get_state(state = None):
    if(state == None):
      if(random.uniform(0.0, 100.0) > 50.0):
        return True
      return False
    else:
      if(random.uniform(0.0, 100.0) > 50.0):
        if(state):
          return False
        return True
      return state
      
      
  def mutate(params, chance=0.01):
    mutant = pf_cycle.random_p(0.0, g.fit['mutate_scale'])
    for i in range(len(params)):
      if(random.uniform(0.0, 1.0) <= chance):
        params[i] = mutant[i]
    return params   
    
    
  def make_fresh():    
    w = g.fit['fresh_ws']
    w_inc = (g.fit['fresh_we'] - g.fit['fresh_ws']) / (g.fit['fresh_size'] // 2 - 1)
    for p in range(g.fit['fresh_size']):
      loop = True
      while(loop):
        try:
          rn = min(len(g.top_parameters) - 1, numpy.floor((len(g.top_parameters) + 1) * random.uniform(0.0, 1.0)**3))
          pbest = numpy.copy(g.top_parameters[rn][1])

          g.pfdata['params']['fresh'][p, :] = pf_parameters.random_p(pbest)
          pf_potential.update(g.pfdata['params']['fresh'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['fresh_rss'][p] = rss
        except:
          pass 

          
        
  # Used to merge parents, fresh and all children into ordered array      
  def merge():
    ps = g.fit['pop_size']
    fs = g.fit['fresh_size']
    cs = g.fit['pop_size'] + 2 * g.fit['fresh_size']
    
  
    # Reset array
    g.pfdata['params']['pop_merged'][:,:] = 0.0
    g.pfdata['params']['pop_merged_rss'][:] = 0.0
    g.pfdata['params']['pop_merged'][0:ps,:] = g.pfdata['params']['pop'][:,:]   
    g.pfdata['params']['pop_merged_rss'][0:ps] = g.pfdata['params']['pop_rss'][:] 
    g.pfdata['params']['pop_merged'][ps:ps+fs,:] = g.pfdata['params']['fresh'][:,:]   
    g.pfdata['params']['pop_merged_rss'][ps:ps+fs] = g.pfdata['params']['fresh_rss'][:] 
    g.pfdata['params']['pop_merged'][ps+fs:ps+fs+cs,:] = g.pfdata['params']['children'][:,:]  
    g.pfdata['params']['pop_merged_rss'][ps+fs:ps+fs+cs] = g.pfdata['params']['children_rss'][:] 
        
    # Sort
    pf_generation.sort_merged()
    
  def sort_merged():
    sort.sort_1d_dp_asc(g.pfdata['params']['pop_merged_rss'][:])
    g.pfdata['params']['pop_merged_rss'] = sort.apply_keytable_1d_dp(g.pfdata['params']['pop_merged_rss'])
    g.pfdata['params']['pop_merged'] = sort.apply_keytable_2d_dp(g.pfdata['params']['pop_merged'])

  
  




















###############################################################################################
#
# Delete some day
#
###############################################################################################    

  def make_fresh_old():    
    w = g.fit['fresh_ws']
    w_inc = (g.fit['fresh_we'] - g.fit['fresh_ws']) / (g.fit['fresh_size'] // 2 - 1)
    for p in range(g.fit['fresh_size']):
      loop = True
      while(loop):   
        if(p % 2 == 0):
          g.pfdata['params']['fresh'][p, :] = pf_cycle.random_p(g.pfdata['top']['p'][0,:], w)       
        else:  
          r = random.randint(0, 9)
          g.pfdata['params']['fresh'][p, :] = pf_cycle.random_p(g.pfdata['top']['p'][r,:], w)              
        loop = False
        try:
          pf_potential.update(g.pfdata['params']['fresh'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['fresh_rss'][p] = rss
        except:
          pass 
      if(p % 2 == 1):
        w = w + w_inc
  