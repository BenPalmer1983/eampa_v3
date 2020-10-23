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
from pot_fit_generation import pf_generation
from pot_fit_extinction import pf_extinction
from pot_fit_enhance import pf_enhance
from pot_fit_parameters import pf_parameters
import matplotlib.pyplot as plt
import time
import random
import hashlib


class pf:
  
  def run():    
    
    # Set up EFS and BP modules
    pf.set_up()
    
    # Make data structure
    g.pfdata = pf_data.make()
        
    # Print start screen
    pf.startup()
    
    # Make Pool
    pf_parameters.make_pool()
   
    # Run
    pf.run_fit()
    
    # Save Results
    pf.save()
    
    
  ######################################################     
  # SET UP - initialise EFS, BP and any other modules
  #          and run first calc
  ###################################################### 
  def set_up():  
    # If a spline fit, convert into a spline with the set number of nodes
    for fn in range(len(g.pot_functions['functions'])):       
      if(g.pot_functions['functions'][fn]['fit_type'] == 1): 
        potential.vary_tabulated_points(fn)
  
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
    # Setup BP    
    bp_calc.init()
    potential.bp_add_potentials()
    b_props.bp_add()
    bp_calc.get_known()
    
    # Rescale density function
    if(g.fit['rescale_density'] == 2):
      rescale_density.run()      
    
    
    
  def startup():
    display.clear() 
    
    
    print("###################################################")      
    print("Starting RSS: ", g.pfdata['rss']['start'])   
    print("###################################################") 
    print("Potential Parameters")
    print("###################################################")   
    potential.print_parameters() 
    
    time.sleep(2.5)

    
  def run_fit():
  

    for c in range(g.fit['cycles']):
      pf_cycle.run()
  
    
     
     
     
     
     
     
     
    
  def get_rss(top=True):  
    rss = rss_calc.run_calc()  
    
    
    g.pfdata['rss']['current'] = rss
    g.pfdata['rss']['counter'] += 1
    
    if(rss is not None):    
      # Hash of Parameters
      h = pf.param_hash(g.pfdata['params']['current'])    
      
      # Store Details
      g.pfdata['rss']['counter_successful'] += 1
      g.pfdata['rss']['since_improvement'] += 1    
      if(g.pfdata['rss']['start'] == None):
        g.pfdata['rss']['start'] = rss
      if(g.pfdata['rss']['best'] == None or rss < g.pfdata['rss']['best']):
        g.pfdata['params']['best'][:] = copy.deepcopy(g.pfdata['params']['current'])
        g.pfdata['best_hash'] = h
        g.pfdata['rss']['best'] = rss
        g.pfdata['bp_best'] = copy.deepcopy(g.bp_results)        
        g.pfdata['rss']['since_improvement'] = 0
        
      # Fill in top 
      if(top):
        pf.top(g.pfdata['params']['current'], rss)

    # DISPLAY
    display.output()    
    return rss
    
    
  def top(p, rss):  
    g.pfdata['top']['counter'] += 1
    if(g.pfdata['top']['filled']):
      if(rss < g.pfdata['top']['rss'][-1]):
        g.pfdata['top']['rss'][-1] = copy.deepcopy(rss)
        g.pfdata['top']['p'][-1,:] = copy.deepcopy(p)
    else:
      for n in range(g.pfdata['top']['size']):
        breakout = False
        if(g.pfdata['top']['rss'][n] == -1.0):
          g.pfdata['top']['rss'][n] = copy.deepcopy(rss)
          g.pfdata['top']['p'][n,:] = copy.deepcopy(p)
          breakout = True
        if(breakout):
          break
      if(breakout == False):
        g.pfdata['top']['filled'] = True
    
    sort.sort_1d_dp_asc(g.pfdata['top']['rss'][:])
    g.pfdata['top']['rss'] = sort.apply_keytable_1d_dp(g.pfdata['top']['rss'])
    g.pfdata['top']['p'] = sort.apply_keytable_2d_dp(g.pfdata['top']['p'])

    #top = {'size': top_size, 'rss': numpy.zeros((top_size,),), 'p': numpy.zeros((top_size, width,),),}
    #pass
  
    """
    p = g.pfdata['params']['current']
    for i in range(10):
      if(rss == g.pfdata['top_ten'][i][0]):
        return 0
      if(g.pfdata['top_ten'][i][0] == None):
        g.pfdata['top_ten'][i][0] = rss
        g.pfdata['top_ten'][i][1] = p
        return 0
      if(rss < g.pfdata['top_ten'][i][0]):
        if(i < 9):
          g.pfdata['top_ten'][i+1:9] = copy.deepcopy(g.pfdata['top_ten'][i:8])
        g.pfdata['top_ten'][i][0] = rss
        g.pfdata['top_ten'][i][1] = p
        return 0
    """
    
    
  def print_top_ten():  
    print()  
    print("=====================")  
    for i in range(10):
      print(i, g.pfdata['top_ten'][i][0])
    print("=====================")  
    print()  
    print()  
    
    
  def param_hash(p):
    hstr = ''
    for i in range(len(p)):
      hstr = hstr + str(p[i])
    return hashlib.md5(hstr.encode()).hexdigest()
    
    
    
    
    
  def save():  
   
    # Load best
    pf_potential.update(g.pfdata['params']['best'][:])
    
    # Display
    g.pfdata['stage'] = 'Finished'
    display.finish()
    
    # Output
    potential_output.full()
    
  
  
  
  
  
  
   
  
  
  
  
  
  
  
  