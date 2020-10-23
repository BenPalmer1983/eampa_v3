######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import matplotlib.pyplot as plt
import time

class efs_calc:

  def run_energy():
  
    print("Calc Energy") 
  
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy()
    efs_calc.output_energy()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  
  def run_energy_force():
  
    print("Calc Energy and Forces") 
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
  
  def run_energy_force_stress():
  
    print("Calc Energy, Forces and Stress") 
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force_stress() 
    efs_calc.output_energy()
    efs_calc.output_forces()
    efs_calc.output_stress()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  

    #efs.max_density_calc()

  def output_energy():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_energies.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_energy[n,:], f_pad)
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  
  
  def output_forces():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_forces.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('FORCE RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('\n')
    for n in range(efs.cc): 
      fh.write('##################\n')
      fh.write('Config ' + str(n) + '\n')
      fh.write('##################\n')
      a = efs.key[n, 0] - 1
      b = efs.key[n, 1]
      for l in range(a, b):
        std.write_file_line(fh, str(efs.labels[l]) + ':', t_pad, efs.config_forces[l,:], f_pad)
      fh.write('\n')
     
      
    fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()



  def output_stress():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['results'] + '/' + 'config_stresses.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('STRESS RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    for n in range(efs.cc):    
      a = n * 3
      std.write_file_line(fh, 'Config ' + str(n+1) + ':', t_pad, efs.config_stresses[a,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+1,:], f_pad)
      std.write_file_line(fh, '', t_pad, efs.config_stresses[a+2,:], f_pad)
      fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  

  def output_rss():
    fh = open(g.dirs['results'] + '/' + 'rss_configs.txt', 'w')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################\n')
    for n in range(efs.cc):
      fh.write(str(efs.config_energy[n,0]) + ' ')
      fh.write(str(efs.config_energy[n,1]) + ' ')
      fh.write(str(efs.config_energy[n,2]) + ' ')
      fh.write(str(efs.config_energy[n,3]) + ' ')
      fh.write(str(efs.config_energy[n,4]) + ' ')
      fh.write(str(efs.config_energy[n,5]) + ' ')
      fh.write('\n')
    fh.write('###############################################################\n')
    fh.close()






  def set_weights():
    # READ INPUT DATA
    efs.set_weights(g.rss_weights['config'], g.rss_weights['energy'], g.rss_weights['force'], g.rss_weights['stress'])
  
  
  
  
  def get_results():
    efs_calc.get_known()    
    # SET UP DICTIONARY
    g.efs_results = {}    
    try:
      g.efs_results['ok'] = True
      g.efs_results['cc'] = int(efs.cc)
      g.efs_results['efs_calculations'] = []
      for cn in range(int(efs.cc)):
        g.efs_results['efs_calculations'].append({'e': None,})
        g.efs_results['efs_calculations'][cn]['e'] = float(efs.config_energy[cn, 2])
    except:
      g.efs_results['ok'] = False
      g.efs_results['cc'] = 0 
      
       
  def get_known():  
    # SET UP DICTIONARY
    g.efs_known = {} 
    try:        
      g.efs_known['ok'] = True
      g.efs_known['cc'] = int(efs.cc)
      g.efs_known['efs_calculations'] = []
      for cn in range(int(efs.cc)):
        g.efs_known['efs_calculations'].append({'e': None,})
        g.efs_known['efs_calculations'][cn]['e'] = float(efs.energies[cn])
    except:
      g.efs_known['ok'] = False
      g.efs_known['cc'] = 0       
       
      
      
  def get_rss():
    try:
      # Make Dictionary
      g.rss['efs'] = {}
      g.rss['efs']['ok'] = True
      g.rss['efs']['cc'] = int(efs.cc)
      g.rss['efs']['energy_rss'] = float(efs.energy_rss)
      g.rss['efs']['force_rss'] = float(efs.force_rss)
      g.rss['efs']['stress_rss'] = float(efs.stress_rss)
      g.rss['efs']['total_rss'] = float(efs.total_rss)
      g.rss['efs']['energy_rss_weighted'] = float(efs.energy_rss_weighted)
      g.rss['efs']['force_rss_weighted'] = float(efs.force_rss_weighted)
      g.rss['efs']['stress_rss_weighted'] = float(efs.stress_rss_weighted)
      g.rss['efs']['total_rss_weighted'] = float(efs.total_rss_weighted)
    except:
      # Make Dictionary
      g.rss['efs'] = {}
      g.rss['efs']['ok'] = False
      g.rss['efs']['cc'] = 0
  
 
    """
    try:
      # Make Dictionary
      g.rss_efs = {}
      g.rss_efs['ok'] = True
      g.rss_efs['cc'] = int(efs.cc)
      g.rss_efs['energy_rss'] = float(efs.energy_rss)
      g.rss_efs['force_rss'] = float(efs.force_rss)
      g.rss_efs['stress_rss'] = float(efs.stress_rss)
      g.rss_efs['total_rss'] = float(efs.total_rss)
      g.rss_efs['energy_rss_weighted'] = float(efs.energy_rss_weighted)
      g.rss_efs['force_rss_weighted'] = float(efs.force_rss_weighted)
      g.rss_efs['stress_rss_weighted'] = float(efs.stress_rss_weighted)
      g.rss_efs['total_rss_weighted'] = float(efs.total_rss_weighted)
    except:
      # Make Dictionary
      g.rss_efs = {}
      g.rss_efs['ok'] = False
      g.rss_efs['cc'] = 0
    """ 
  
