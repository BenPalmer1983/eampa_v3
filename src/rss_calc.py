######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from potential import potential
import matplotlib.pyplot as plt
import time


class rss_calc:

  def run():
    print("Calc RSS") 
    
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    efs_calc.set_weights()
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
    
    # Setup BP
    bp_calc.init()
    bp_calc.set_weights()
    potential.bp_add_potentials()
    b_props.bp_add()

    # Run EFS and BP
    rss = rss_calc.run_calc()
    

    # Output to File
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    
    # Plots
    #potential.plot_fortran_potentials()
    #potential.plot_python_potentials()
    potential.plot_python_potentials(g.dirs['wd'] + '/plots/potential_python')
    potential.plot_fortran_potentials(g.dirs['wd'] + '/plots/potential_fortran')
    
    b_props.bp_output()
    b_props.bp_eos_plot()
    
    
    
    print('')     
    print('CONFIGS')    
    for n in range(efs.cc):    
      print('Config ' + str(n+1) + ':', efs.config_energy[n,2], efs.energies[n], (efs.config_energy[n,2]-efs.energies[n])**2)
    print('All configs:                   ' + str(efs.total_rss))
    print('All configs (energy):          ' + str(efs.energy_rss))
    print('All configs (force):           ' + str(efs.force_rss))
    print('All configs (stress):          ' + str(efs.stress_rss))
    print('All configs weighted:          ' + str(efs.total_rss_weighted))
    print('All configs weighted (energy): ' + str(efs.energy_rss_weighted))
    print('All configs weighted (force):  ' + str(efs.force_rss_weighted))
    print('All configs weighted (stress): ' + str(efs.stress_rss_weighted))
    

    
    print('')   
    for bp_id in range(bp.bp_configs_count):  
      print('BP') 
      print('alat:', bp.calc_alat[bp_id], bp.known_alat[bp_id], (bp.calc_alat[bp_id] - bp.known_alat[bp_id])**2)
      print('v0:', bp.calc_v0[bp_id])
      print('e0:', bp.calc_e0[bp_id], bp.known_e0[bp_id], (bp.calc_e0[bp_id] - bp.known_e0[bp_id])**2)
      print('b0:', bp.calc_b0[bp_id], bp.known_b0[bp_id], (bp.calc_b0[bp_id] - bp.known_b0[bp_id])**2)
      print("Calculated Stiffness Matrix (GPA)")
      for i in range(6):
        print(160.230732254e0 * bp.calc_ec[bp_id,i,:])
      print("Known Stiffness Matrix (GPA)")
      for i in range(6):
        print(160.230732254e0 * bp.known_ec[bp_id,i,:])
  
    
    print('')
    print('RSS: ' + str(rss))
    
    
    print('')
    print(g.rss)
    print('')
    
    
    
    
    
  # ASSUMES EFS AND BP ALREADY SET UP
  def run_calc(): 
    # START TIME
    s = time.time()  
  
    # EFS CONFIG COUNT
    try:
      efs_cc = int(efs.cc)
    except:  
      efs_cc = 0
      
    # BP CONFIG COUNT
    try:
      bp_cc = int(bp.cc)
    except:  
      bp_cc = 0
  
    # RUN CALCULATIONS
    if(g.rss_weights['config'] != 0.0 and efs_cc > 0):
      efs.rss_calc()
      efs_calc.get_rss()
    if(g.rss_weights['bp'] != 0.0 and bp_cc > 0):
      bp.energy()
      bp.calculate_bp()  
  
    # LOAD CALCULATED RESULTS
    efs_calc.get_results()
    bp_calc.get_results()
  
    # LOAD RSS RESULTS
    efs_calc.get_rss()
    bp_calc.get_rss()

    # SUM WEIGHTED RSS
    rss = 0.0
    if(g.rss['efs']['ok']):
      rss = rss + g.rss['efs']['total_rss_weighted']
    if(g.rss['bp']['ok']):
      rss = rss + g.rss['bp']['total_rss_weighted'] 
    
    # Increment counter
    if(g.rss['counter'] == None):
      g.rss['counter'] = 1
    else:
      g.rss['counter'] = g.rss['counter'] + 1
    
    # If bad result, return None
    if(numpy.isnan(rss)):
      # RETURN
      g.rss['current'] = None
      g.rss['log'].append(None)
      g.rss['counter'] = g.rss['counter'] + 1
      return None
      
    # LOG
    #rss = {'current': None, 'best': None, 'counter': None, 'since_improvement': None, 'log': []} 
    g.rss['current'] = rss
    g.rss['log'].append(rss)
    
    # KEEP LOG OF BEST RSS
    if(g.rss['since_improvement'] == None):
      g.rss['since_improvement'] = 0
    g.rss['since_improvement'] = g.rss['since_improvement'] + 1
    if(g.rss['best'] == None or g.rss['current'] < g.rss['best']):
      g.rss['best'] = g.rss['current']
      g.rss['since_improvement'] = 0
      
      g.efs_results_best = copy.deepcopy(g.efs_results)
      g.bp_results_best = copy.deepcopy(g.bp_results)
      
    # END TIME
    e = time.time()  
    
    dt = e - s
    
    
    #efs_cc bp_cc
    g.benchmark['total_atoms'] = g.benchmark['total_atoms'] + (bp.total_atoms + efs.total_atoms)
    g.benchmark['total_interactions'] = g.benchmark['total_interactions'] + (bp.l_nl_size + efs.l_nl_size)
    g.benchmark['configs'] =  g.benchmark['configs'] + efs_cc + bp_cc
    g.benchmark['total_time'] = g.benchmark['total_time'] + dt
    g.benchmark['atomspersec'] =  g.benchmark['total_atoms'] / g.benchmark['total_time'] 
    g.benchmark['interationspersec'] =  g.benchmark['total_interactions'] / g.benchmark['total_time'] 
    g.benchmark['configspersec'] =  g.benchmark['configs'] / g.benchmark['total_time'] 
    
    # RETURN
    return rss
    
    
    

    
  def reset_rss(): 
    g.rss = {'current': None, 'best': None, 'counter': None, 'since_improvement': None, 'log': [], 'efs': {'ok': False, 'cc': 0,}, 'bp': {'ok': False, 'cc': 0,}}
    
    
    
    