######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy

"""
1. init()
This prepares the arrays and clear any data

2. add configs for BP testing (e.g. BCC, FCC, and so on)
This makes the configs, ghost config and the neighbour list

3. add potentials
Update with the latest potential

4. calculate energy (no forces, no stresses)
Calculates and saves the energy for each


"""



class bp_calc:

  def run():  
  
    print("Calc Bulk Properties") 
    
    # Setup BP
    bp_calc.init()
    
    # Load potentials    
    potential.bp_add_potentials()
    
    # Add required BP structures (add, make ghost configs, make NLs)
    b_props.bp_add()
        
    # Calculate energies of configurations for all structures
    bp.energy() 
   
        
    # Calculate BP
    bp.calculate_bp()   
    
    
    #print(bp.rss)
    
    
    # Output to Terminal
    b_props.bp_output_terminal()
    
    # Output to File
    #b_props.bp_output()
    
    # Plots
    #b_props.bp_eos_plot()
    #potential.plot_fortran_potentials()
    #potential.plot_python_potentials()
    
    
    
    
  def init():
    # Log
    g.log_fh.write('BP Allocated Memory\n')
    try:
      mem = str(g.inp['mem']['bp'])
    except:
      mem = "500MB"
    g.log_fh.write(mem + '\n')
    mem = std.mem_value(mem)    
    g.memory['bp']['c'] = int(1 * (mem / 7840))
    g.memory['bp']['g'] = int(12 * (mem / 7840))
    g.memory['bp']['nl'] = int(100 * (mem / 7840))
    g.log_fh.write(str(g.memory['bp']['c']) + " " + str(g.memory['bp']['g']) + " " + str(g.memory['bp']['nl']) + '\n')
    
    # Initialise
    bp.init(g.memory['bp']['c'], g.memory['bp']['g'], g.memory['bp']['nl'])
    
    
    
    
  def set_weights():
    bp.set_rss_total(g.rss_weights['bp'])
    bp.set_rss_a0(g.rss_weights['a0'])
    bp.set_rss_e0(g.rss_weights['e0'])
    bp.set_rss_b0(g.rss_weights['b0'])
    bp.set_rss_ec(g.rss_weights['ec'])
    bp.set_rss_g(g.rss_weights['g'])
    bp.set_rss_e(g.rss_weights['e'])
    bp.set_rss_v(g.rss_weights['v'])
    
        
        
  def get_results(): 
    bp_calc.get_known()  
    # Make Dictionary
    g.bp_results = {}    
    try:             
      g.bp_results['ok'] = True
      g.bp_results['cc'] = int(bp.cc)
      g.bp_results['bp_calculations'] = []
      g.bp_results['input'] = {}
      bp_id = 0
      while(bp.bp_keys_i[bp_id,0] > -1):
        g.bp_results['input'][bp_id] = g.bp_ids[bp_id+1]
        g.bp_results['bp_calculations'].append({'a0': None, 'e0': None, 'b0': None, 'ec': None, 'g': None, 'e': None, 'v': None,'b0_gpa': None,'ec_gpa': None,})
        if(bp.known_set[bp_id, 0] == 1):
          g.bp_results['bp_calculations'][bp_id]['a0'] = float(bp.calc_alat[bp_id])          
        if(bp.known_set[bp_id, 1] == 1):
          g.bp_results['bp_calculations'][bp_id]['e0'] = float(bp.calc_e0[bp_id])        
        if(bp.known_set[bp_id, 2] == 1):
          g.bp_results['bp_calculations'][bp_id]['b0'] = float(bp.calc_b0[bp_id])   
          g.bp_results['bp_calculations'][bp_id]['b0_gpa'] = 160.230732254 * float(bp.calc_b0[bp_id])       
        if(bp.known_set[bp_id, 3] == 1):
          g.bp_results['bp_calculations'][bp_id]['ec'] = numpy.zeros((6,6,),)
          g.bp_results['bp_calculations'][bp_id]['ec_gpa'] = numpy.zeros((6,6,),)
          for i in range(6):
            for j in range(6):
              g.bp_results['bp_calculations'][bp_id]['ec'][i,j] = float(bp.calc_ec[bp_id, i, j])
              g.bp_results['bp_calculations'][bp_id]['ec_gpa'][i,j] = 160.230732254 * float(bp.calc_ec[bp_id, i, j])
        if(bp.known_set[bp_id, 4] == 1):
          g.bp_results['bp_calculations'][bp_id]['g'] = float(bp.calc_g[bp_id])
        if(bp.known_set[bp_id, 5] == 1):
          g.bp_results['bp_calculations'][bp_id]['e'] = float(bp.calc_e[bp_id])
        if(bp.known_set[bp_id, 6] == 1):
          g.bp_results['bp_calculations'][bp_id]['v'] = float(bp.calc_v[bp_id])
        bp_id = bp_id + 1        
    except:
      g.bp_results['ok'] = False
      g.bp_results['cc'] = 0 


  def get_known(): 
    # Make Dictionary
    g.bp_known = {}   
    try:    
      g.bp_known['ok'] = True
      g.bp_known['cc'] = int(bp.cc)
      g.bp_known['bp_calculations'] = {}
      g.bp_known['input'] = {}
      bp_id = 0
      while(bp.bp_keys_i[bp_id,0] > -1):   
        g.bp_known['input'][bp_id] = g.bp_ids[bp_id+1]
        g.bp_known['bp_calculations'][bp_id] = {'a0': None, 'e0': None, 'b0': None, 'ec': None, 'g': None, 'e': None, 'v': None,'b0_gpa': None,'ec_gpa': None,}
        
        
        if(bp.known_set[bp_id, 0] == 1):
          g.bp_known['bp_calculations'][bp_id]['a0'] = float(bp.known_alat[bp_id])
        if(bp.known_set[bp_id, 1] == 1):
          g.bp_known['bp_calculations'][bp_id]['e0'] = float(bp.known_e0[bp_id])
        if(bp.known_set[bp_id, 2] == 1):
          g.bp_known['bp_calculations'][bp_id]['b0'] = float(bp.known_b0[bp_id])
          g.bp_known['bp_calculations'][bp_id]['b0_gpa'] = 160.230732254 * float(bp.known_b0[bp_id])
        if(bp.known_set[bp_id, 3] == 1):
          g.bp_known['bp_calculations'][bp_id]['ec'] = numpy.zeros((6,6,),)
          g.bp_known['bp_calculations'][bp_id]['ec_gpa'] = numpy.zeros((6,6,),)
          for i in range(6):
            for j in range(6):
              g.bp_known['bp_calculations'][bp_id]['ec'][i,j] = float(bp.known_ec[bp_id, i, j])
              g.bp_known['bp_calculations'][bp_id]['ec_gpa'][i,j] = 160.230732254 * float(bp.known_ec[bp_id, i, j])
        if(bp.known_set[bp_id, 4] == 1):
          g.bp_known['bp_calculations'][bp_id]['g'] = float(bp.known_g[bp_id])
        if(bp.known_set[bp_id, 5] == 1):
          g.bp_known['bp_calculations'][bp_id]['e'] = float(bp.known_e[bp_id])
        if(bp.known_set[bp_id, 6] == 1):
          g.bp_known['bp_calculations'][bp_id]['v'] = float(bp.known_v[bp_id])
        bp_id = bp_id + 1  
    except:
      g.bp_known['ok'] = False
      g.bp_known['cc'] = 0 
    
    
    
  def get_rss():  
     
    try:
      # Make Dictionary
      g.rss['bp'] = {}
      g.rss['bp']['ok'] = True  
      g.rss['bp']['cc'] = int(bp.cc)
      # Totals
      g.rss['bp']['total_rss'] = float(bp.rss_total_rss)
      g.rss['bp']['total_rss_weighted'] = float(bp.rss_total_rss_w)
      # Individual RSS
      g.rss['bp']['a0'] = float(bp.rss_by_type[0])
      g.rss['bp']['e0'] = float(bp.rss_by_type[1])
      g.rss['bp']['b0'] = float(bp.rss_by_type[2])
      g.rss['bp']['ec'] = float(bp.rss_by_type[3])
      g.rss['bp']['g'] = float(bp.rss_by_type[4])
      g.rss['bp']['e'] = float(bp.rss_by_type[5])    
      g.rss['bp']['v'] = float(bp.rss_by_type[6])
      # Weighted RSS
      g.rss['bp']['a0_weighted'] = float(bp.rss_by_type_w[0])
      g.rss['bp']['e0_weighted'] = float(bp.rss_by_type_w[1])
      g.rss['bp']['b0_weighted'] = float(bp.rss_by_type_w[2])
      g.rss['bp']['ec_weighted'] = float(bp.rss_by_type_w[3])
      g.rss['bp']['g_weighted'] = float(bp.rss_by_type_w[4])
      g.rss['bp']['e_weighted'] = float(bp.rss_by_type_w[5])    
      g.rss['bp']['v_weighted'] = float(bp.rss_by_type_w[6])
    except:
      # Make Dictionary
      g.rss['bp'] = {}
      g.rss['bp']['ok'] = False
      g.rss['bp']['cc'] = 0 
           
      
      
      
      
      
    """   
    try:
      # Make Dictionary
      g.rss_bp = {}
      g.rss_bp['ok'] = True  
      g.rss_bp['cc'] = int(bp.cc)
      # Totals
      g.rss_bp['total_rss'] = float(bp.rss_total_rss)
      g.rss_bp['total_rss_weighted'] = float(bp.rss_total_rss_w)
      # Individual RSS
      g.rss_bp['a0'] = float(bp.rss_by_type[0])
      g.rss_bp['e0'] = float(bp.rss_by_type[1])
      g.rss_bp['b0'] = float(bp.rss_by_type[2])
      g.rss_bp['ec'] = float(bp.rss_by_type[3])
      g.rss_bp['g'] = float(bp.rss_by_type[4])
      g.rss_bp['e'] = float(bp.rss_by_type[5])    
      g.rss_bp['v'] = float(bp.rss_by_type[6])
      # Weighted RSS
      g.rss_bp['a0_weighted'] = float(bp.rss_by_type_w[0])
      g.rss_bp['e0_weighted'] = float(bp.rss_by_type_w[1])
      g.rss_bp['b0_weighted'] = float(bp.rss_by_type_w[2])
      g.rss_bp['ec_weighted'] = float(bp.rss_by_type_w[3])
      g.rss_bp['g_weighted'] = float(bp.rss_by_type_w[4])
      g.rss_bp['e_weighted'] = float(bp.rss_by_type_w[5])    
      g.rss_bp['v_weighted'] = float(bp.rss_by_type_w[6])
    except:
      # Make Dictionary
      g.rss_bp = {}
      g.rss_bp['ok'] = False
      g.rss_bp['cc'] = 0 
    """  
    
