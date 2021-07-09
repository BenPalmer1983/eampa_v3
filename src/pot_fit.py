######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from potential import potential
from gd import gd
from rescale_density import rescale_density
from pot_fit_start import pf_start
from pot_fit_steps import pf_steps
from pot_fit_genetic import pf_genetic
from pot_fit_data import pf_data
from pot_fit_cycle import pf_cycle
from pot_fit_potential import pf_potential
from pot_fit_generation import pf_generation
from pot_fit_extinction import pf_extinction
from pot_fit_enhance import pf_enhance
from pot_fit_parameters import pf_parameters
from pot_fit_setop import pf_setop
from pot_fit_gd import pf_gd
from pot_fit_ng import pf_ng
from pot_fit_sa import pf_sa
from pot_fit_neldermead import pf_nm
from pot_fit_bh import pf_bh
from pot_fit_cg import pf_cg
from pot_fit_bfgs import pf_bfgs
from pot_fit_gsa import pf_gsa
from pot_fit_saa import pf_saa
from pot_fit_sps import pf_sps
from pot_fit_saved import pf_saved
from pot_fit_random import pf_saved
from pot_fit_top import pf_top
from pot_fit_pgradient import pf_top
from pot_fit_display import pf_display
from pot_fit_display_simple import pf_display_simple
from pot_fit_init import pf_init
from pot_fit_save import pf_save
from pot_fit_final import pf_final
import matplotlib.pyplot as plt
import time
import random
import hashlib


class pf:

  start_time = 0.0
  end_time = 0.0  
  rss_plot_data = []

  def run():   
  
    # Start Time
    pf.start_time = time.time()

    pf.pbest_dir = g.dirs['wd'] + '/fitting/saved_potentials/0000_pbest/'
    std.make_dir(pf.pbest_dir)
  
    # Log
    main.log_title("Potential Fit")
      
    # Set up EFS and BP modules and g.pfdata
    pf_init.run()

    # Create Plot Dirs
    std.make_dir(g.dirs['wd'] + '/fitting')
    pf.fh = open(g.dirs['wd'] + '/fitting/fitting_summary.txt', 'w')


    # Run once with initial parameters
    g.pfdata['stage'] = 'START - INITIAL PARAMETERS'
    g.pfdata['stage_brief'] = 'START'
    rss = pf.get_rss()


    for fit in g.fit:
      pf.fit = fit
      if(pf.fit['type'].lower() == "random"):        
        pf_random.run(fit)
      elif(pf.fit['type'].lower() == "sa"): 
        pf_sa.run(fit)
      elif(pf.fit['type'].lower() == "gsa"): 
        pf_gsa.run(fit)
      elif(pf.fit['type'].lower() == "ga"): 
        pf_genetic.run(fit)
      elif(pf.fit['type'].lower() == "ng"): 
        pf_ng.run(fit)
      elif(pf.fit['type'].lower() == "sps"): 
        pf_sps.run(fit)
      elif(pf.fit['type'].lower() == "nm"): 
        pf_neldermead.run(fit)
      elif(pf.fit['type'].lower() == "bfgs"): 
        pf_bfgs.run(fit)
      elif(pf.fit['type'].lower() == "bh"): 
        pf_bh.run(fit)
      elif(pf.fit['type'].lower() == "cg"): 
        pf_cg.run(fit)


    pf_save.top("BEST_PARAMETERS")

   
    potential.update(g.pfdata['p']['best'] )
    rss = pf.get_rss()

    pf_top.save(g.dirs['wd'] + '/save', 'top_parameters.txt')

    # Output
    #potential_output.full()
          
    # End Time
    pf.end_time = time.time()
  
    # Log fit time
    main.log("Fit time: ", str(pf.end_time - pf.start_time))

    pf_save.rss_plot_full()



  
  
  
  
  
  ######################################################     
  # SET UP - initialise EFS, BP and any other modules
  #          and run first calc
  ###################################################### 

    
  def get_rss(save_in_top=True, quiet=False):  
    rss = rss_calc.run_calc(save_in_top)  
        
    g.pfdata['rss']['current'] = rss
    g.pfdata['rss']['counter'] += 1

    # If successful
    if(rss is not None and not quiet):   
    
      # Store Details
      g.pfdata['bp']['current'] = g.bp_results.copy()
      g.pfdata['rss']['counter_successful'] += 1
      g.pfdata['rss']['since_improvement'] += 1    
      if(g.pfdata['rss']['start'] == None):
        g.pfdata['rss']['start'] = rss
      if(g.pfdata['rss']['best'] == None or rss < g.pfdata['rss']['best']):
        g.pfdata['p']['best'] = numpy.copy(potential.pot['p'])
        g.pfdata['rss']['best'] = rss
        g.pfdata['bp']['best'] = g.bp_results.copy()    
        g.pfdata['efs']['best'] = g.efs_results.copy()     
        g.pfdata['rss']['since_improvement'] = 0       

        potential.save(pf.pbest_dir, g.pfdata['p']['best']) 
        pf_save.save_summary(pf.pbest_dir, g.pfdata['p']['best'], g.pfdata['rss']['best'], g.pfdata['bp']['best'], g.pfdata['bp']['known'],  g.pfdata['efs']['best'], g.pfdata['efs']['known'])
        potential.plot(pf.pbest_dir)

      # Save in top
      if(save_in_top):
        list_len = g.fitting['top_parameters']

        # COPY ARRAYS
        p_now = numpy.copy(potential.pot['p'])
        efs_results = g.efs_results.copy()
        bp_results = g.bp_results.copy()
        rss_details = g.rss.copy()

        if(len(g.top_parameters) == 0):
          g.top_parameters.append([rss, p_now, efs_results, bp_results, bp.max_density, efs.max_density, rss_details])
        else:    
          i = 0
          while(i<len(g.top_parameters)):
            if(rss < g.top_parameters[i][0]):
              break
            elif(rss == g.top_parameters[i][0] and g.top_parameters[i][1].all() == p_now.all()):
              i = list_len
              break
            else:
              i = i + 1
          if(i < list_len):
            g.top_parameters.insert(i, [rss, p_now, efs_results,bp_results, bp.max_density, efs.max_density, rss_details])
          if(len(g.top_parameters) > list_len):
            g.top_parameters.pop()    



    # DISPLAY
    if(not quiet):
      pf_display.output()    
    #for i in range(min(len(g.top_parameters),10)):
    #  print(g.top_parameters[i][0])
    return rss
    
    
    
    

  def summary_line(opt_type, t_start, t_end, rss_start, rss_end):
    
    t_taken = str('{:6.3f}'.format(t_end - t_start))
    while(len(t_taken)<16):
      t_taken = t_taken + " "
    
    rss_start = str('{:12.4e}'.format(rss_start))
    while(len(rss_start)<16):
      rss_start = rss_start + " "
    
    rss_end = str('{:12.4e}'.format(rss_end))
    while(len(rss_end)<16):
      rss_end = rss_end + " "

    opt_type = str(opt_type)
    while(len(opt_type)<30):
      opt_type = opt_type + " "

    pf.fh.write(opt_type + rss_start + "-->  " + rss_end + "  " + t_taken + "\n")
    
    
  def save_rss_plot_data(t_start, name, data):
    plot_data = numpy.asarray(data)
    pf.rss_plot_data.append([t_start, name, data])



