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


class pf_genetic:

  start_time = 0

  def run():    
    if(g.fit['cycles'] == 0 or g.fit['gens'] == 0):
      return 0
    
  
    pf_genetic.start_time = time.time()
  
    main.log_title("Genetic Fit")    
    
    # Print start screen
    #pf_genetic.startup()
    
    # Output
    potential.plot_python_potentials(g.dirs['wd'] + '/plots/pots/start_potential')
   
    # Run
    pf_genetic.run_fit()
    
    # Save Results
    pf_genetic.save()
    
    

    
    
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
  
    # Loop through cycles
    for c in range(g.fit['cycles']):
      pf_cycle.run()
  
    
    
  def save():  
   
    main.log_hr()
    main.log("End of fit")
    main.log("Time spent in potfit:   " + str(time.time() - pf_genetic.start_time))
    main.log("RSS Counter:            " + str(g.rss['counter']))
    main.log("Configs:                " + str(g.benchmark['configs']))
    main.log("Total atoms:            " + str(g.benchmark['total_atoms']))
    main.log("Total Interactions:     " + str(g.benchmark['total_interactions']))
    main.log("Time:                   " + str(g.benchmark['total_time']))
    main.log("Configs/Sec:            " + str(g.benchmark['configspersec']))
    main.log("Atoms/Sec:              " + str(g.benchmark['atomspersec']))
    main.log("Interactions/Sec:       " + str(g.benchmark['interationspersec']))
    
    main.log_hr()
    main.log("Parameters")    
    for i in range(len(g.pfdata['params']['best'])):
      main.log("P" + str(i) + " " + str(g.pfdata['params']['best'][i]))
    main.log_hr()
   
   
   
    # Load best
    pf_potential.update(g.pfdata['params']['best'][:])
    

  
  
  
  
  
   
  
  
  
  
  
  
  
  