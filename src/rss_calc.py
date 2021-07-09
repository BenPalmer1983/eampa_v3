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

    # Make dir
    rss_calc.output_dir()   

    
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    efs_calc.set_weights()
    potential.load_to_efs()
    configs.efs_add_config()           # Add configs

    # Setup BP
    bp_calc.init()
    bp_calc.set_weights()
    potential.load_to_bp()
    b_props.bp_add()


    # Run EFS and BP
    rss = rss_calc.run_calc()
    
    
    efs_calc.output()
    efs_calc.save_to_file(g.dirs['wd'] + '/rss')


    b_props.bp_output_terminal()


    print('')     
    print('CONFIGS')   
    print('             e known         e calc         e_rss         e_rss_w       f_rss       f_rss_w       s_rss       s_rss_w')      
    for n in range(efs.cc):    
      print('Config ' + str(n+1) + ':', end="")
      print('{:18.6f}'.format(efs.energies[n]), end="")
      print('{:18.6f}'.format(efs.config_energy[n,2]), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,0])), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,4])), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,1])), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,5])), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,2])), end="")
      print('{:18.6f}'.format(float(efs.config_rss[n,6])), end="")
      print('')     

    e_rss = sum(efs.config_rss[:,0])
    f_rss = sum(efs.config_rss[:,1])
    s_rss = sum(efs.config_rss[:,2])
    t_rss = e_rss + f_rss + s_rss

    e_rss_w = sum(efs.config_rss[:,4])
    f_rss_w = sum(efs.config_rss[:,5])
    s_rss_w = sum(efs.config_rss[:,6])
    t_rss_w = e_rss_w + f_rss_w + s_rss_w

    print()
    print('All configs:                   ' + str(t_rss))
    print('All configs (energy):          ' + str(e_rss))
    print('All configs (force):           ' + str(f_rss))
    print('All configs (stress):          ' + str(s_rss))
    print('All configs weighted:          ' + str(t_rss_w))
    print('All configs weighted (energy): ' + str(e_rss_w))
    print('All configs weighted (force):  ' + str(f_rss_w))
    print('All configs weighted (stress): ' + str(s_rss_w))

   
    print('')   
    for bp_id in range(bp.bp_configs_count):  
      print('')   
 
      print('========================================================================================') 
      print('========================================================================================') 
      print('BP  ' + str(bp_id)) 
      print('========================================================================================') 
      print('========================================================================================') 

      rss_calc.print_line('a0', bp.known_alat[bp_id], bp.calc_alat[bp_id])
      rss_calc.print_line('v0', None, bp.calc_v0[bp_id])
      rss_calc.print_line('e0', bp.known_e0[bp_id], bp.calc_e0[bp_id])
      rss_calc.print_line('b0', bp.known_b0[bp_id], bp.calc_b0[bp_id])
      print('')   

      print("Calculated Stiffness Matrix (GPA)")
      for i in range(6):
        for j in range(6):
          print(str('{:14.6f}'.format(float(160.230732254e0 * bp.calc_ec[bp_id,i,j]))), end="")
        print()
      print("Known Stiffness Matrix (GPA)")
      for i in range(6):
        for j in range(6):
          print(str('{:14.6f}'.format(float(160.230732254e0 * bp.known_ec[bp_id,i,j]))), end="")
        print()
        #print(160.230732254e0 * bp.known_ec[bp_id,i,:])

      print('')   
      print('Other calculated properties') 
      print('===========================')     
      print('{:30s}'.format('Bulk Modulus B0 (Reuss):'), end=" ")
      print('{:14.6f}'.format(bp.calc_b0_r[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_b0_r_gpa[bp_id]),") GPA", sep="")
      print('{:30s}'.format('Bulk Modulus B0 (Voigt):'), end=" ")
      print('{:14.6f}'.format(bp.calc_b0_v[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_b0_v_gpa[bp_id]),") GPA", sep="")

      print('{:30s}'.format('Tetragonal Shear:'), end=" ")
      print('{:14.6f}'.format(bp.calc_cubic_tetragonal_shear[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_cubic_tetragonal_shear[bp_id]),") GPA", sep="")

      print('{:30s}'.format('Shear Modulus G:'), end=" ")
      print('{:14.6f}'.format(bp.calc_cubic_shear_modulus[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_cubic_shear_modulus_gpa[bp_id]),") GPA", sep="")
      print('{:30s}'.format('Shear Modulus G (Reuss):'), end=" ")
      print('{:14.6f}'.format(bp.calc_g_r[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_g_r_gpa[bp_id]),") GPA", sep="")
      print('{:30s}'.format('Shear Modulus G (Voigt):'), end=" ")
      print('{:14.6f}'.format(bp.calc_g_v[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_g_v_gpa[bp_id]),") GPA", sep="")

      print('{:30s}'.format('Young\'s Modulus E:'), end=" ")
      print('{:14.6f}'.format(bp.calc_e[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_e_gpa[bp_id]),") GPA", sep="")

      print('{:30s}'.format('Poisson Ratio'), end=" ")
      print('{:14.6f}'.format(bp.calc_v[bp_id]))

      print('{:30s}'.format('Cachy Pressure:'), end=" ")
      print('{:14.6f}'.format(bp.calc_cubic_cauchy_pressure[bp_id]), end="    ")
      print("(", '{:14.6f}'.format(bp.calc_cubic_cauchy_pressure_gpa[bp_id]),") GPA", sep="")

      print('{:30s}'.format('Melting Point:'), end=" ")
      print('{:14.6f}'.format(bp.calc_melting[bp_id]), "K", sep="")

      print('{:30s}'.format('Cubic Stability:'), end=" ")
      print('{:14.6f}'.format(bp.calc_cubic_stability[bp_id]), " (1=Stable, 0=Unstable)", sep="")
      print('{:30s}'.format('Stability:'), end=" ")
      print('{:14.6f}'.format(bp.calc_stability[bp_id]), " (1=Stable, 0=Unstable)", sep="")
      print('')         
      print('')   
  
    print('RSS:')    
    print('a0:', g.rss['bp']['a0'], "   w: ",g.rss['bp']['a0_weighted'])
    print('e0:', g.rss['bp']['e0'], "   w: ",g.rss['bp']['e0_weighted'])
    print('b0:', g.rss['bp']['b0'], "   w: ",g.rss['bp']['b0_weighted'])
    print('ec:', g.rss['bp']['ec'], "   w: ",g.rss['bp']['ec_weighted'])
    print('g:', g.rss['bp']['g'], "   w: ",g.rss['bp']['g_weighted'])
    print('e:', g.rss['bp']['e'], "   w: ",g.rss['bp']['e_weighted'])
    print('v:', g.rss['bp']['v'], "   w: ",g.rss['bp']['v_weighted'])
    print('')
    print('')
    print('')
    print('')
    print('RSS: ' + str(rss))
    print('')
    print('')
    print('Max Density:', bp.max_density, efs.max_density)
    print('')
    print('')
    #print(g.rss)
    print('')
    
    potential.plot(g.dirs['wd'] + '/rss/potential_function_plots')

    std.make_dir(g.dirs['wd'] + '/rss/eos_ec')  
    b_props.bp_eos_plot(g.dirs['wd'] + '/rss/eos_ec')

    # Plots
    #potential.plot_fortran_potentials()
    #potential.plot_python_potentials()
    #potential.plot_python_potentials(g.dirs['wd'] + '/rss/plots/potential_python')
    #potential.plot_fortran_potentials(g.dirs['wd'] + '/rss/plots/potential_fortran')
    #b_props.bp_output()
    #b_props.bp_eos_plot(g.dirs['wd'] + '/rss/plots')


    
    

    

    


  def output_dir():
    std.make_dir(g.dirs['wd'] + '/rss')  
    std.make_dir(g.dirs['wd'] + '/rss/plots')
    
    
  # ASSUMES EFS AND BP ALREADY SET UP
  def run_calc(save_in_top=True, log=False): 
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
    efs.rss_calc()
    bp.energy()
    bp.calculate_bp()  

    # GET RESULTS
    efs_calc.get_results()
    efs_calc.get_rss()   
    bp_calc.get_results()
    bp_calc.get_rss()

    # Log
    if(log == True):
      rss_calc.log()
 
    # Save residual
    g.rss['residual'] = numpy.asarray(g.rss['residual'])
  

    # SUM WEIGHTED RSS
    rss = 0.0
    if(g.rss['efs']['ok']):
      rss = rss + g.rss['efs']['total_rss_weighted']
    if(g.rss['bp']['ok']):
      rss = rss + g.rss['bp']['total_rss_weighted'] 

    #print("Python EFS total_rss_weighted", g.rss['efs']['total_rss_weighted'])
    #print("Python BP total_rss_weighted", g.rss['bp']['total_rss_weighted'])
    #print( g.rss['efs']['total_rss_weighted'] + g.rss['bp']['total_rss_weighted'])
    #print(sum(g.rss['residual']))  

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
      
    
    # KEEP LOG OF BEST RSS
    if(g.rss['since_improvement'] == None):
      g.rss['since_improvement'] = 0
    g.rss['since_improvement'] = g.rss['since_improvement'] + 1
    g.rss['current_bp_max_density'] = numpy.copy(bp.max_density)
    g.rss['current_efs_max_density'] = numpy.copy(efs.max_density)
    if(g.rss['best'] == None or g.rss['current'] < g.rss['best']):
      main.log('Best rss: ' + str(g.rss['best']) + ' (since last: ' + str(g.rss['since_improvement']) + ')')
      g.rss['best_bp_max_density'] = numpy.copy(bp.max_density)
      g.rss['best_efs_max_density'] = numpy.copy(efs.max_density)
      g.rss['best'] = g.rss['current']
      g.rss['since_improvement'] = 0
      
      g.efs_results_best = numpy.copy(g.efs_results)
      g.bp_results_best = numpy.copy(g.bp_results)

    # END TIME
    e = time.time()      
    dt = e - s

    
    #efs_cc bp_cc
    g.benchmark['dt'] = dt
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
    g.rss = {'current': None, 
             'best': None, 
             'counter': None, 
             'since_improvement': None, 
             'top': [], 
             'efs': {'ok': False, 'cc': 0,}, 
             'bp': {'ok': False, 'cc': 0,}
            }


  def print_line(label, known, calc):
    while(len(label)<14):
      label = label + " "

    if(known == None):
      print(label, end="")
      print("              ", end="")
      print(str('{:14.6f}'.format(calc)), end="")
      print()
    elif(calc == None):
      print(label, end="")
      print(str('{:14.6f}'.format(known)), end="")
      print("              ", end="")
      print()
    elif(not(known == None and calc == None)):
      print(label, end="")
      print(str('{:14.6f}'.format(known)), end="")
      print(str('{:14.6f}'.format(calc)), end="")
      print(str('{:14.6f}'.format((calc - known)**2)), end="")
      print()



  def log():
    std.make_dir(g.dirs['wd'] + '/logs') 
    fh.open(g.dirs['wd'] + '/logs/efs_calc.txt', 'a')  
    fh.close()
    print(g.efs_results )
    