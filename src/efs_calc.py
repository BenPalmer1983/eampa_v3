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
    efs_calc.output_dir()
  
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy()
    efs_calc.output()
    efs_calc.save_to_file()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  
  def run_energy_force():
  
    print("Calc Energy and Forces") 
    efs_calc.output_dir()
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force() 
    efs_calc.output()
    efs_calc.save_to_file()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
  
  
  def run_energy_force_stress():
  
    print("Calc Energy, Forces and Stress") 
    efs_calc.output_dir()
    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    efs.energy_force_stress() 
    configs.efs_results()              # Load Results
    efs_calc.output()
    efs_calc.save_to_file()
    #efs_calc.output_rss()
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()

    #efs.max_density_calc()

  def output_dir():
    std.make_dir(g.dirs['wd'] + '/efs')


  def output():
  
    t_pad = 12
    f_pad = 18
    margin = 30
    halfwidth = 30

    print()

    for n in range(efs.cc):    

      # Config calcs      
      nat = efs.key[n,1] - efs.key[n,0] + 1
      a = efs.key[n, 0]
      b = efs.key[n, 1]

      e_on = efs.key[n, 10]
      f_on = efs.key[n, 11]
      s_on = efs.key[n, 12]
      #print(efs.key[n,:])

      print("Config " + str(n+1) + '    ' + g.configs['configs'][n]['file_path'])
      print('###############################################################')
      print(std.pad("Atom Count:", margin) + str(nat))

      print(std.pad("Energy (known/calculated):", margin) + 
            std.pad(efs.energies[n], halfwidth) + 
            std.pad(efs.config_energy[n,2] , halfwidth))

      print(std.pad("Pair:", margin) + 
            std.pad("", halfwidth) + 
            std.pad(efs.config_energy[n,0] , halfwidth))

      print(std.pad("Embedding:", margin) + 
            std.pad("", halfwidth) + 
            std.pad(efs.config_energy[n,1] , halfwidth))


      print()

  def save_to_file(out_dir=None):

    if(out_dir == None):
      out_dir = g.dirs['wd'] + '/efs'
  
    t_pad = 12
    f_pad = 18
    margin = 20
    halfwidth = 50
  
    fh = open(out_dir + '/efs_results.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('FULL RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('\n')    
    for n in range(efs.cc):    

      # Config calcs      
      nat = efs.key[n,1] - efs.key[n,0] + 1
      a = efs.key[n, 0]
      b = efs.key[n, 1]

      e_on = efs.key[n, 10]
      f_on = efs.key[n, 11]
      s_on = efs.key[n, 12]
      #print(efs.key[n,:])

      fh.write("Config " + str(n+1) + '    ')
      fh.write(g.configs['configs'][n]['file_path'] + '\n')
      fh.write('###############################################################\n')
      fh.write('\n')

      

      fh.write(std.pad("Atom Count:", margin))
      fh.write(str(nat))
      fh.write('\n')

      fh.write(std.pad("", margin))
      fh.write(std.pad("KNOWN", halfwidth))
      fh.write(std.pad("CALCULATED", halfwidth))
      fh.write('\n')
      
      # Energy
      fh.write(std.pad("ENERGY:", margin))
      fh.write('\n')

      fh.write(std.pad("(total)", margin))
      fh.write(std.pad(efs.energies[n], halfwidth))
      fh.write(std.pad(efs.config_energy[n,2] , halfwidth))
      fh.write('\n')

      fh.write(std.pad("(pair)", margin))
      fh.write(std.pad("", halfwidth))
      fh.write(std.pad(efs.config_energy[n,0], halfwidth))
      fh.write('\n')

      fh.write(std.pad("(embedding)", margin))
      fh.write(std.pad("", halfwidth))
      fh.write(std.pad(efs.config_energy[n,1], halfwidth))
      fh.write('\n')

      fh.write(std.pad("EPA:", margin))
      fh.write('\n')

      fh.write(std.pad("(total)", margin))
      fh.write(std.pad(efs.energies[n] / nat, halfwidth))
      fh.write(std.pad(efs.config_energy[n,5], halfwidth))
      fh.write('\n')

      fh.write(std.pad("(pair)", margin))
      fh.write(std.pad("", halfwidth))
      fh.write(std.pad(efs.config_energy[n,3], halfwidth))
      fh.write('\n')

      fh.write(std.pad("(embedding)", margin))
      fh.write(std.pad("", halfwidth))
      fh.write(std.pad(efs.config_energy[n,4], halfwidth))
      fh.write('\n')
      fh.write('\n')

      

      fh.write(std.pad("Atom Coords:", margin))
      fh.write('\n')
      nn = 0
      for cn in range(a, b, 1):
        fh.write(std.pad(str(nn) + ":", 8))
        fh.write(std.pad(labels.get(efs.labels[cn]) + " [" + str(efs.labels[cn]) + "]", 20))
        fh.write(std.pad('{:6.3f}'.format(efs.coords[cn,0]), 14))
        fh.write(std.pad('{:6.3f}'.format(efs.coords[cn,1]), 14))
        fh.write(std.pad('{:6.3f}'.format(efs.coords[cn,2]), 14))

        if(f_on == 1):
          fh.write("  #  ")
          fh.write(std.pad('{:12.3e}'.format(float(efs.forces[n,nn,0])), 14))
          fh.write(std.pad('{:12.3e}'.format(float(efs.forces[n,nn,1])), 14))
          fh.write(std.pad('{:12.3e}'.format(float(efs.forces[n,nn,2])), 14))
          fh.write("  #  ")
          fh.write(std.pad('{:12.3f}'.format(float(efs.config_forces[n,nn,0])), 14))
          fh.write(std.pad('{:12.3f}'.format(float(efs.config_forces[n,nn,1])), 14))
          fh.write(std.pad('{:12.3f}'.format(float(efs.config_forces[n,nn,2])), 14))
          

        fh.write('\n')
        nn = nn + 1
      fh.write('\n')


      fh.write('\n')     
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()





  def output_energy():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['wd'] + '/efs/' + 'config_energies.txt', 'w')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('ENERGY RESULTS\n')
    fh.write('Config Count: ' + str(efs.cc) + '\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.write('\n')    
    for n in range(efs.cc):    
      fh.write("Config " + str(n+1) + '    ')
      fh.write(g.configs['configs'][n]['file_path'] + '\n')
      std.write_file_line(fh, "", t_pad, efs.config_energy[n,:], f_pad)
      fh.write('\n')     
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()
  
  
  def output_forces():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['wd'] + '/efs/' + 'config_forces.txt', 'w')
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
      for l in range(0, efs.key[n, 19]):
        std.write_file_line(fh, str(efs.labels[l]) + ':', t_pad, efs.config_forces[n, l,:], f_pad)
      fh.write('\n')
     
      
    fh.write('\n')
    fh.write('###############################################################')
    fh.write('###############################################################\n')
    fh.close()



  def output_stress():
  
    t_pad = 12
    f_pad = 18
  
    fh = open(g.dirs['wd'] + '/efs/' + 'config_stresses.txt', 'w')
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
    fh = open(g.dirs['wd'] + '/efs/' + 'rss_configs.txt', 'w')
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

    # Bias so worse RSS if the density is out of range
    """
    f = 1.0  
    if(bp.max_density < g.rss_max_density['min']):
      f = 1.0 + (100.0 * (g.rss_max_density['min'] - bp.max_density))**4
    if(bp.max_density > g.rss_max_density['max']):
      f = 1.0 + (g.rss_max_density['scale_factor'] * (bp.max_density - 0.8))**g.rss_max_density['scale_exponent']
    """
    f = 1.0  
    if(bp.max_density == 0.0):
      f = g.rss_max_density['zero_density_factor']  
    
    g.rss['residual'] = []

    try:
      # Make Dictionary
      g.rss['efs'] = {}
      g.rss['efs']['ok'] = True
      g.rss['efs']['cc'] = int(efs.cc)
      g.rss['efs']['energy_rss'] = float(efs.energy_rss)
      g.rss['efs']['force_rss'] = float(efs.force_rss)
      g.rss['efs']['stress_rss'] = float(efs.stress_rss)
      g.rss['efs']['total_rss'] = float(efs.total_rss)
      g.rss['efs']['energy_rss_weighted'] = f * float(efs.energy_rss_weighted)
      g.rss['efs']['force_rss_weighted'] = f * float(efs.force_rss_weighted)
      g.rss['efs']['stress_rss_weighted'] = f * float(efs.stress_rss_weighted)
      g.rss['efs']['total_rss_weighted'] = f * float(efs.total_rss_weighted)
    except:
      # Make Dictionary
      g.rss['efs'] = {}
      g.rss['efs']['ok'] = False
      g.rss['efs']['cc'] = 0
  

    for i in range(efs.residuals_size):
      g.rss['residual'].append(f * float(efs.residuals[i]))


    
