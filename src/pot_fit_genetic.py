######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from potential import potential
from display import display
from gd import gd
from rescale_density import rescale_density
from pot_fit_data import pf_data
from pot_fit_cycle import pf_cycle
from pot_fit_potential import pf_potential
from pot_fit_extinction import pf_extinction
from pot_fit_enhance import pf_enhance
from pot_fit_parameters import pf_parameters
import matplotlib.pyplot as plt
import time
import random
import hashlib


class pf_genetic:

  count = 0
  start_time = 0
  variation_factor = 1.0

  def run(fit):    
    if(fit['gens'] == 0):
      return 0
    
    #
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best']  
    pf_genetic.start_time = time.time()
    pf_genetic.count = pf_genetic.count + 1
    pf_genetic.rss_plot = []
    pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, start_best_rss])
    pf_genetic.rss_best = start_best_rss

    g.pfdata['stage'] = 'Genetic Algorithm ' + str(pf_genetic.count)
    g.pfdata['stage_brief'] = 'GA' + str(pf_genetic.count)
  
    main.log_title("Genetic Fit")   
    
    # Load from input
    pf_genetic.width = g.pfdata['psize']
    pf_genetic.pop_size = fit['pop_size']
    pf_genetic.fresh_size = fit['fresh_size']
    pf_genetic.generations = fit['gens']
    pf_genetic.no_clone_var = fit['no_clone_var']
    pf_genetic.gen_variation_multiplier = fit['gen_variation_multiplier']


 
    # Force to be even
    if(pf_genetic.pop_size % 2 != 0):
      pf_genetic.pop_size = pf_genetic.pop_size + 1
    if(pf_genetic.fresh_size % 2 != 0):
      pf_genetic.fresh_size = pf_genetic.fresh_size + 1

    # Calculate other sizes
    pf_genetic.children_size = pf_genetic.pop_size + 2 * pf_genetic.fresh_size
    pf_genetic.merge_size = 2 * pf_genetic.pop_size + 3 * pf_genetic.fresh_size
    
    # Set Up Pop + Fresh arrays
    pf_genetic.pop = numpy.zeros((pf_genetic.pop_size, pf_genetic.width,),)
    pf_genetic.pop_rss = numpy.zeros((pf_genetic.pop_size,),)
    pf_genetic.fresh = numpy.zeros((pf_genetic.fresh_size, pf_genetic.width,),)
    pf_genetic.fresh_rss = numpy.zeros((pf_genetic.fresh_size,),)
    pf_genetic.children = numpy.zeros((pf_genetic.children_size, pf_genetic.width,),)
    pf_genetic.children_rss = numpy.zeros((pf_genetic.children_size,),)
    pf_genetic.pop_merged = numpy.zeros((pf_genetic.merge_size, pf_genetic.width,),)
    pf_genetic.pop_merged_rss = numpy.zeros((pf_genetic.merge_size,),)
  
    # Parents array
    pf_genetic.parents = numpy.arange(pf_genetic.pop_size)

    # Initialise population
    g.pfdata['stage'] = 'Genetic Fit ' + str(pf_genetic.count) + ' Initialising'
  
    # Copy out parameters from "top_parameters"
    start_parameters = []
    if(len(g.top_parameters) > 0):
      for i in range(len(g.top_parameters)):
        new_p = numpy.copy(g.top_parameters[i][1][:])
        start_parameters.append(new_p)
    else:
      p = potential.get_start_parameters()
      new_p = numpy.copy(p)
      start_parameters.append(new_p)

    # Load start parameters into pop, and fill remaining with random parameters
    n = 0
    for pn in range(pf_genetic.pop_size):
      loop = True
      while(loop): 
        n = n + 1
        if(n <= len(start_parameters)):
          pf_genetic.pop[pn, :] = start_parameters[n-1][:]
        else:
          pf_genetic.pop[pn, :] = potential.random(g.pfdata['p']['best'], 1.0, False)        
        # Try - if it fails or rss == None, try next
        potential.update(pf_genetic.pop[pn, :])
        try:
          # Update
          rss = pf.get_rss()
          if(rss is not None):
            loop = False
            pf_genetic.pop_rss[pn] = rss
        except:
          pass
        if(rss < pf_genetic.rss_best):
          pf_genetic.rss_best = rss
          pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, pf_genetic.rss_best])

    ###################################
    # LOOP THROUGH GENERATIONS
    ###################################

    pf_genetic.generation = 0
    for gen in range(pf_genetic.generations):
      pf_genetic.generation = pf_genetic.generation + 1
      pf_genetic.run_gen()

    pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, pf_genetic.rss_best])

    # End
    pf.summary_line("GENETIC ALGORITHM " + str(pf_genetic.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("GENETIC_ALGORITHM_" + str(pf_genetic.count))
    pf_save.rss_plot("GENETIC_ALGORITHM_" + str(pf_genetic.count), pf_genetic.rss_plot)
    pf.save_rss_plot_data(t_start, "GENETIC_ALGORITHM_" + str(pf_genetic.count), pf_genetic.rss_plot)
  


  def run_gen():

    ss = 'Genetic Fit ' + str(pf_genetic.count) + ' Gen ' + str(pf_genetic.generation) + ' '  

    # Shuffle parents array
    numpy.random.shuffle(pf_genetic.parents) 

    c = 0
    # Breed Population
    g.pfdata['stage'] = ss + 'Breed Population'
    for p in range(pf_genetic.pop_size // 2):  
      loop = True
      while(loop):
        loop = pf_genetic.breed_event(p, c, 'p+p')      
      c = c + 2

    # Make Fresh  
    g.pfdata['stage'] = ss + 'Initialising Fresh Population'   
    pf_genetic.make_fresh()


    # Shuffle parents array
    numpy.random.shuffle(pf_genetic.parents) 

    # Breed Population-Fresh 
    g.pfdata['stage'] = ss + 'Breed With Fresh Population'  
    for p in range(pf_genetic.fresh_size):  
      loop = True
      while(loop):
        loop = pf_genetic.breed_event(p, c, 'p+f')      
      c = c + 2

    # Merge populations and select top to form next generation   
    g.pfdata['stage'] = ss + 'Merge Populations'   
    pf_genetic.merge()


    # Load best back into population
    pf_genetic.pop[:, :] = pf_genetic.pop_merged[0:pf_genetic.pop_size, :]
    pf_genetic.pop_rss[:] = pf_genetic.pop_merged_rss[0:pf_genetic.pop_size]
   
    # Change variation for next generation
    pf_genetic.variation_factor = pf_genetic.variation_factor * pf_genetic.gen_variation_multiplier




  def breed_event(p, c, opt):
    pc = pf_genetic.width
    
    pa = pf_genetic.parents[p]
    if(opt == 'p+p'):  
      pb = pf_genetic.parents[p + pf_genetic.pop_size // 2]
      pb_array = 'pop'
    if(opt == 'p+f'): 
      pb = p
      pb_array = 'fresh'
    
    
    ca = c
    cb = c + 1
    
    # Breed 
    state = pf_genetic.get_state()
    for i in range(pc):
      state = pf_genetic.get_state(state)
      if(state):
        pf_genetic.children[ca, i] = pf_genetic.pop[pa, i]
        if(pb_array == 'pop'):
          pf_genetic.children[cb, i] = pf_genetic.pop[pb, i]
        elif(pb_array == 'fresh'):
          pf_genetic.children[cb, i] = pf_genetic.fresh[pb, i]
      else:     
        if(pb_array == 'pop'):
          pf_genetic.children[ca, i] = pf_genetic.pop[pb, i]
        elif(pb_array == 'fresh'):
          pf_genetic.children[ca, i] = pf_genetic.fresh[pb, i]
        pf_genetic.children[cb, i] = pf_genetic.pop[pa, i]
      
    # Mutate
    #g.pfdata['params']['children'][ca, :] = pf_genetic.mutate(g.pfdata['params']['children'][ca, :], g.fit['mutate_chance'])
    #g.pfdata['params']['children'][cb, :] = pf_genetic.mutate(g.pfdata['params']['children'][cb, :], g.fit['mutate_chance'])
    

    # No clones - Child A
    for pn in range(pf_genetic.pop_size):
      if((pf_genetic.children[ca, :].all() == pf_genetic.pop[pn,:]).all()):
        r = numpy.random.rand(pc)   
        pf_genetic.children[ca, :] = pf_parameters(pf_genetic.children[ca, :],  pf_genetic.no_clone_var)
        break
        
    # No clones - Child B
    for pn in range(pf_genetic.pop_size):
      if((pf_genetic.children[cb, :].all() == pf_genetic.pop[pn,:]).all()):
        r = numpy.random.rand(pc)   
        pf_genetic.children[cb, :] = pf_parameters(pf_genetic.children[cb, :],  pf_genetic.no_clone_var)
        break

    # Check the children actually give valid parameters
    loop = False     
    potential.update(pf_genetic.children[ca, :])
    rss_a = pf.get_rss()
    if(rss_a is None):
      loop = True 
      
    potential.update(pf_genetic.children[cb, :])
    rss_b = pf.get_rss()
    if(rss_b is None):
      loop = True
      
     
    if(loop == False):
      pf_genetic.children_rss[ca] = rss_a
      pf_genetic.children_rss[cb] = rss_b
      
    if(rss_a != None and rss_b != None):
      if(rss_a < pf_genetic.rss_best):
        pf_genetic.rss_best = rss_a
        pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, pf_genetic.rss_best])
      if(rss_b < pf_genetic.rss_best):
        pf_genetic.rss_best = rss_b
        pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, pf_genetic.rss_best])


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



  def make_fresh():    
    for pn in range(pf_genetic.fresh_size):
      loop = True
      pbest = numpy.copy(g.top_parameters[0][1])
      while(loop):
        try:
          rn = int(min(len(g.top_parameters) - 1, numpy.floor((len(g.top_parameters) + 1) * random.uniform(0.0, 1.0)**3)))
          p_trial = numpy.copy(g.top_parameters[rn][1])
          p_new = potential.random(p_trial, pf_genetic.variation_factor, False)
          rss = pf.get_rss()
          if(rss is not None):
            loop = False
        except:
          pass 
      pf_genetic.fresh[pn, :] = p_new
      pf_genetic.fresh_rss[pn] = rss

      if(rss < pf_genetic.rss_best):
        pf_genetic.rss_best = rss
        pf_genetic.rss_plot.append([time.time() - pf_genetic.start_time, pf_genetic.rss_best])


  # Used to merge parents, fresh and all children into ordered array      
  def merge():
    ps = pf_genetic.pop_size
    fs = pf_genetic.fresh_size
    cs = pf_genetic.pop_size + 2 * pf_genetic.fresh_size
      
    # Reset array
    pf_genetic.pop_merged[:,:] = 0.0
    pf_genetic.pop_merged[0:ps,:] = pf_genetic.pop[:,:]   
    pf_genetic.pop_merged[ps:ps+fs,:] = pf_genetic.fresh[:,:]   
    pf_genetic.pop_merged[ps+fs:ps+fs+cs,:] = pf_genetic.children[:,:]  
    pf_genetic.pop_merged_rss[:] = 0.0
    pf_genetic.pop_merged_rss[0:ps] = pf_genetic.pop_rss[:] 
    pf_genetic.pop_merged_rss[ps:ps+fs] = pf_genetic.fresh_rss[:] 
    pf_genetic.pop_merged_rss[ps+fs:ps+fs+cs] = pf_genetic.children_rss[:] 
        
    # Sort
    pf_genetic.sort_merged()
    
  def sort_merged():
    sort.sort_1d_dp_asc(pf_genetic.pop_merged_rss[:])
    pf_genetic.pop_merged_rss = sort.apply_keytable_1d_dp(pf_genetic.pop_merged_rss)
    pf_genetic.pop_merged = sort.apply_keytable_2d_dp(pf_genetic.pop_merged)

  
  


  