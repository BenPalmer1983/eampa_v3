######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from potential import potential
from display import display
from gd import gd
from rescale_density import rescale_density
import matplotlib.pyplot as plt
import time
import random





"""
    'f_on': 1,  
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=EMBE, 3=DENS
    'f_group': 1,
    'r_cut': 6.5,
    'file': None,
    'function_type': 0,       # 1 tab, 2 analytic
    'f_points': None,         # READ IN TO PYTHON
    'a_type': '',
    'f': None,
    'a_params': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
"""


class pf:

  
  # Array Sizes
  pop_size = 50  
  fresh_size = 10  
  pop_overflow = 10
  pop_size_d = None
  fresh_size_d = None 
  arr_len = None 
  
  # Parameter arrays
  ps = None
  p_var = None
  p_width = 0

  # Rescale Density  0.0 to 1.0 (max is approximated)
  rescale_on = True
  
  # Variance Coeffs
  vc_extinction = 0.2

  
  ncv = 0.05
  tab_nodes = 10
  generations = 200
  spline_generations = 20
  breed_switch = 0.5
  cycles = 1

  extinction_counter = 0
  extinction_percentage = 50
  extinction_frequency = 3
  extinction_threshold_t = 0.0
  extinction_threshold_k = 0.0

  # Record alat, e0, b0
  bulk_properties = None

  # PROGRESS COUNTERS
  stage = 'INITIALISING'
  start_time = 0.0
  time_per_calc = 0.0
  time_estimate = 0.0
  time_remaining = 0.0
  since_improvement = 0
  this_gen = 0
  this_cycle = 0
  total_generations = 0
  last_rss = 0.0
  best_rss = 0.0
  calc_counter = 0
  calc_counter_expected = 0
  progress = 0

  # Potential Output
  output_counter = 0
  
  # Result
  result_ps = None
  result_rss = None



  data = {
  'start_time': 0.0,
  'end_time': None,
  'cycle_start_time': None,
  'cycles': 0,
  'cycle': 0,  
  'spline_cycles': 0,
  'spline_cycle': 0,
  'generations': 0,
  'generation': 0,
  'spline_generations': 0,
  'spline_generation': 0,
  'population_size': 0,
  'population': 0,
  'fresh_size': 0,
  'population_fs': 0,
  'extinction_percentage': 0.0,
  'extinction_frequency': 0.0,
  'enhance_frequency': 0.0,
  'since_improvement': 0,  
  'weight_config': 1.0,
  'weight_force': 1.0,
  'weight_stress': 1.0,
  'weight_energy': 1.0,
  'weight_bp': 1.0,
  'weight_a0': 1.0,
  'weight_b0': 1.0,
  'weight_surface_energy': 1.0,
  'weight_vacancy': 1.0,
  'a0_known': 0.0,
  'a0_calc': 0.0,
  'a0_rss': 0.0,
  'e0_known': 0.0,
  'e0_calc': 0.0,
  'e0_rss': 0.0,
  'config_rss': 0.0,
  'config_rss_weighted': 0.0,
  'rss_counter': 0,
  'best_rss': None,
  'last_rss': None,
  }






  def run():    
    
    # Read parameter variance
    try:
      pf.pv_in = g.inp['fit']['variance']
    except:
      pass
      
    # Read from input into globals
    try:
      g.fit['generations'] = g.inp['fit']['gens']
    except:
      pass
    try:
      g.fit['spline_generations'] = g.inp['fit']['sgens']
    except:
      pass
    try:
      g.fit['spline_cycles'] = g.inp['fit']['scycles']
    except:
      pass
    try:
      g.fit['pop_size'] = g.inp['fit']['ps']
    except:
      pass
    try:
      g.fit['fresh_size'] = g.inp['fit']['fs']
    except:
      pass
    try:
      g.fit['extinction_percentage'] = g.inp['fit']['ep']
    except:
      pass
    try:
      g.fit['extinction_frequency'] = g.inp['fit']['ef']
    except:
      pass
    try:
      g.fit['enhance_frequency'] = g.inp['fit']['ehf']
    except:
      pass
    try:
      g.fit['cycles'] = g.inp['fit']['cycles']
    except:
      pass
      
      
    # Set pf globals
    pf.cycles = g.fit['cycles']
    pf.generations = g.fit['generations']
    pf.spline_cycles = g.fit['spline_cycles']
    pf.spline_generations = g.fit['spline_generations']
    pf.pop_size = g.fit['pop_size']
    pf.fresh_size = g.fit['fresh_size']
    pf.extinction_percentage = g.fit['extinction_percentage']
    pf.extinction_frequency = g.fit['extinction_frequency']
    pf.enhance_frequency = g.fit['enhance_frequency']
    
    pf.data['cycles'] = g.fit['cycles']
    pf.data['spline_cycles'] = g.fit['spline_cycles']
    pf.data['generations'] = g.fit['generations']
    pf.data['spline_generations'] = g.fit['spline_generations']
    pf.data['population_size'] = g.fit['pop_size']
    pf.data['fresh_size'] = g.fit['fresh_size']
    pf.data['extinction_percentage'] = g.fit['extinction_percentage']
    pf.data['extinction_frequency'] = g.fit['extinction_frequency']
    pf.data['enhance_frequency'] = g.fit['enhance_frequency']
   
    
    # Set Up EFS and BP
    pf.set_up()
    
    # Run start up screen
    pf.startup()
    
    # Start fit
    pf.fit()
    
   
  ##################################################     
  # SET UP
  ##################################################  
  def set_up():  
    # If a spline fit, convert into a spline 
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
    
    # Rescale density function
    rescale_density.run()      
      
      
      
      
      
  def startup():
  
    os.system('cls' if os.name == 'nt' else 'clear')    
    start_rss = rss_calc.run_calc()
    
    print("###################################################")      
    print("Starting RSS: ", start_rss)   
    print("###################################################") 
    print("Potential Parameters")
    print("###################################################")   
    pf.print_parameters() 
    
    time.sleep(2.5)
    
    
    
  def print_parameters():  
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("fn: " + str(fn))
        print("type: spline")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("Parameter " + str(i) + ": ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][0,i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][1,i], end='')
          print()
        #pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size'] 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        print("fn: " + str(fn))
        print("type: analytic")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("Parameter " + str(i) + ": ", end='')
          print(g.pot_functions['functions'][fn]['a_params'][i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][0,i], end='')
          print(", ", end='')
          print(g.pot_functions['functions'][fn]['fit_parameters'][1,i], end='')
          print()
      else:
        print("fn: " + str(fn))
        print("type: no variance")

    
   
    
  def get_rss():        
    pf.calc_counter = pf.calc_counter  + 1
    rss = rss_calc.run_calc()    
    return rss
   
    

    
        
  def update_potential(p):    
    # Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE    
        ###NEED TO REDO MAYBE??###
      
        # Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size']          
        # LOAD ORIGINAL
        g.pot_functions['functions'][fn]['points'] = numpy.copy(g.pot_functions['functions'][fn]['points_original'])        
        # VARY SPLINE
        potential.vary_tabulated_points(fn, p[a:b])
        # Update a
        a = b        
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        # Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)
        a = b    
      
    # Rescale density functions 0.0 to 1.0  
    if(pf.rescale_on):
      rescale_density.run()

    
    # Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials



    
    
  def fit():
  
    # Clear Screen
    os.system('cls' if os.name == 'nt' else 'clear') 
    print("START")
  
    # Init counters
    pf.start_time = time.time()
    pf.since_improvement = 0
    pf.improvement_time = time.time()
    pf.this_gen = 0
    pf.extinction_counter = 0
    pf.extinction_threshold_k = 0.0
    pf.extinction_threshold_t = 0.0
    pf.this_cycle = 0
    pf.total_generations = 0
    pf.last_rss = 0.0
    pf.best_rss = 0.0
    pf.calc_counter = 0    
    pf.estimate_calc_counter()
    
    rss_calc.reset_rss()
    
    
    pf.cycle_ps = None  
    c = 1    
    for i in range(pf.cycles):
      pf.data['cycle'] = c
      pf.cycle(c, i, 1)
      c = c + 1
    
    
    exit()
    
    
    
    pf.cycle_ps = None    
    if(pf.spline_cycles>0):
      pf.spline_cycles_prep()
      
    for i in range(pf.spline_cycles):
      pf.spline_cycles_node_var(c)     # Set how much the nodes vary
      pf.cycle(c, i, 2)
      c = c + 1
      
      

    
  def spline_cycles_prep():
    
    # Number of Nodes
    pf.pw_per_f = 20
  
    # Load Best
    pf.update_potential(pf.ps[-5,:])  
    
    # Prep each with a spline
    yvar = numpy.zeros((pf.pw_per_f,),)  
    for fn in range(len(g.pot_functions['functions'])):
      potential.vary_tabulated_points(fn, yvar)
   
    # Update potential fit settings
    for fn in range(len(g.pot_functions['functions'])): 
      pw = 0
      # Replace original points with new points
      g.pot_functions['functions'][fn]['points_original'] = numpy.copy(g.pot_functions['functions'][fn]['points'])  
      # fit type to spline
      g.pot_functions['functions'][fn]['fit_type'] = 1  
      g.pot_functions['functions'][fn]['fit_size'] = pf.pw_per_f    
      g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((2, pf.pw_per_f, ),)
      g.pot_functions['functions'][fn]['fit_mult'] = 1.0
      pw = pw + pf.pw_per_f
      
      
  def spline_cycles_node_var(c):
    pw = 0
    for fn in range(len(g.pot_functions['functions'])):  
      g.pot_functions['functions'][fn]['fit_parameters'][0,:] = -10 * 0.1**c   # Lower
      g.pot_functions['functions'][fn]['fit_parameters'][1,:] = 10 * 0.1**c    # Upper
      g.pot_functions['functions'][fn]['fit_parameters'][:,-1] = 0.0           # Don't move last node
      pw = pw + pf.pw_per_f
      
      
      
      
      
      
      

  def cycle(c, c_id, ct_id):
   
    """
    -1  Initial input parameters
    -2  Mutation
    -3  child a
    -4  child b
    -5  best
    """
    pf.data['cycle_start_time'] = time.time()
    
    # Counter
    pf.this_cycle = pf.this_cycle + 1
    pf.this_gen = 0
       

    # Start Cycle   (calculate width, 
    # Pop size
    pf.pop_size_d = 2 * pf.pop_size
    pf.fresh_size_d = 2 * pf.fresh_size    
    pf.arr_len = pf.pop_size_d + pf.fresh_size_d + pf.pop_overflow
    pf.child_len = pf.pop_size_d + pf.fresh_size_d + pf.fresh_size_d
    
    # Find width  
    pf.p_width = 0
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size'] 
        print(fn, g.pot_functions['functions'][fn]['fit_size'] )
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        pf.p_width = pf.p_width + g.pot_functions['functions'][fn]['fit_size']

    # MAKE ARRAYS
    pf.p_var = numpy.zeros((2, pf.p_width,),)
    pf.ps = numpy.zeros((pf.arr_len, pf.p_width,),)
    pf.rss = numpy.zeros((pf.arr_len,),)   
    pf.p_children = numpy.zeros((pf.child_len, pf.p_width,),)
    pf.rss_children = numpy.zeros((pf.child_len,),)
    pf.rss_list = numpy.zeros((pf.pop_size_d + pf.child_len,),)
    pf.p_temp = numpy.zeros((pf.pop_size_d, pf.p_width,),)
    pf.rss_temp = numpy.zeros((pf.pop_size_d,),)
    pf.p_top = numpy.zeros((int(pf.pop_size_d / 4), pf.p_width,),)
    pf.rss_sort = numpy.zeros((pf.pop_size_d,),)
    
    # SET Variation 
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        pf.p_var[0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        pf.p_var[1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        pf.p_var[0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        pf.p_var[1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  
         
    # SAVE Starting Parameters
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        pf.ps[-1,a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        pf.ps[-1,a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b   
    
    if(pf.cycle_ps is None):
      pf.cycle_ps = numpy.zeros((max(pf.cycles+2, pf.spline_cycles+2), pf.p_width,),)
      pf.cycle_ps[0,:] = pf.ps[-1,:]
        
    # SET -1.0 for best rss   
    pf.rss[-5] = -1.0
    
    
    
    ##############################################
    # CREATE INITIAL POPULATION - FIRST HALF
    ##############################################
    
    pf.data['stage'] = 'Initialising population - First Half'
    for p in range(pf.pop_size):
      if(p == 0):
        pf.ps[p,:] = pf.cycle_ps[c_id,:]
      else:
        # Calculate and save RSS
        pf.ps[p,:] = pf.random_p(0.0, 1.0)
        
      loop = True
      while(loop): 
        # Update potential, run and save rss
        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss() 
        if(not numpy.isnan(pf.rss[p])):
          loop = False
          pf.check_improvement(pf.ps[p,:], pf.rss[p])
        else:
          pf.ps[p,:] = pf.random_p(0.0, 1.0)
      print(pf.ps[p,:])

    # Look in a wider range

    # Create initial population - Second half    
    pf.data['stage'] = 'Initialising population - Second Half'
    m = 0.5
    m_inc = (10.0-0.5) / (pf.pop_size - 1)
    for p in range(pf.pop_size, pf.pop_size_d):
      loop = True
      while(loop):  
        pf.ps[p,:] = pf.random_p(0.0, m)
        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss() 
        if(not numpy.isnan(pf.rss[p])):
          loop = False
          pf.check_improvement(pf.ps[p,:], pf.rss[p])
          m = m + m_inc

    
    ###################################
    # LOOP THROUGH GENERATIONS
    ###################################
    
    gens = pf.generations
    if(c > 1):    
      gens = pf.spline_generations
    pf.since_improvement = 0
    for gen in range(gens):
      pf.data['stage'] = 'Looping population'
      pf.data['generation'] = gen + 1
      
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)

      
      
      #######################
      # Parents + Parents
      #######################
      
      ca = 0
      # Loop through population
      for p in range(pf.pop_size):   
        
        # If nan, keep trying
        loop = True
        retry = 0
        while(loop):  
          pa = parents[p]
          pb = parents[p + pf.pop_size]
        
          # Breed
          pf.breed(pa, pb, ca, ca + 1)
          # Run calculations
          pf.update_potential(pf.p_children[ca, :])        
          pf.rss_children[ca] = pf.get_rss()  
          if(not numpy.isnan(pf.rss_children[ca])):
            pf.check_improvement(pf.p_children[ca,:], pf.rss_children[ca])      
          
          pf.update_potential(pf.p_children[ca+1, :])        
          pf.rss_children[ca+1] = pf.get_rss()  
          if(not numpy.isnan(pf.rss_children[ca+1])):
            pf.check_improvement(pf.p_children[ca+1,:], pf.rss_children[ca+1])            
          if(not numpy.isnan(pf.rss_children[ca]) and not numpy.isnan(pf.rss_children[ca+1])):  
            loop = False
          retry = retry + 1
          if(retry > 4):
            pf.ps[pa,:] = pf.random_p(0.0, m)
            pf.ps[pb,:] = pf.random_p(0.0, m)
            
            
                
        # Increment
        ca = ca + 2
 
 
      
      #######################
      # Parents + Fresh
      #######################
      
      # Pick parents
      parents = numpy.arange(pf.pop_size_d)
      numpy.random.shuffle(parents)
      
      pf.data['stage'] = 'Looping population + fresh pool'
      # Breed random parents with fresh pool
      pa_n = 0
      m_inc = (2.0-0.5) / (pf.fresh_size - 1)
      for p_fresh in range(pf.pop_size_d, pf.pop_size_d + pf.fresh_size_d):
        # Pick Parent
        pa = parents[pa_n]
        pa_n = pa_n + 1
        # If nan, keep trying
        loop = True
        while(loop): 

          rn = numpy.random.randint(pf.fresh_size) + 1
          m = rn * m_inc 
        
          # Make Fresh
          pf.ps[p_fresh,:] = pf.random_p(pf.ps[-5,:], m)
        
          # Breed
          pf.breed(pa, p_fresh, ca, ca + 1)
        
          # Run calculations
          pf.update_potential(pf.p_children[ca, :])        
          pf.rss_children[ca] = pf.get_rss()   
          if(not numpy.isnan(pf.rss_children[ca])):
            pf.check_improvement(pf.p_children[ca, :], pf.rss_children[ca])  
          
          pf.update_potential(pf.p_children[ca+1, :])        
          pf.rss_children[ca + 1] = pf.get_rss()  
          if(not numpy.isnan(pf.rss_children[ca+1])):
            pf.check_improvement(pf.p_children[ca+1, :], pf.rss_children[ca+1])
                   
          if(not numpy.isnan(pf.rss_children[ca]) and not numpy.isnan(pf.rss_children[ca+1])):  
            loop = False

         
        # Increment
        ca = ca + 2
        
        
      
      
      #######################
      # Select Best
      #######################
        
      pf.rss_list[:pf.pop_size_d] = pf.rss[:pf.pop_size_d]
      pf.rss_list[pf.pop_size_d:] = pf.rss_children[:]

      rss_sorted = sort.sort_1d_dp_asc(pf.rss_list)
      rt = rss_sorted[pf.pop_size_d-1]


      # SELECT BEST FROM PARENTS AND CHILDREN FOR NEXT GEN
      n = 0
      for i in range(pf.pop_size_d):
        if(pf.rss[i] < rt):
          pf.p_temp[n,:] = pf.ps[i,:]
          pf.rss_temp[n] = pf.rss[i]
          n = n + 1          
      for i in range(len(pf.rss_children)):
        if(pf.rss_children[i] < rt):
          pf.p_temp[n,:] = pf.p_children[i,:]
          pf.rss_temp[n] = pf.rss_children[i]
          n = n + 1
      for i in range(pf.pop_size_d):
        if(pf.rss[i] == rt and n < pf.pop_size_d):
          pf.p_temp[n,:] = pf.ps[i,:]
          pf.rss_temp[n] = pf.rss[i]
          n = n + 1      
      for i in range(len(pf.rss_children)):
        if(pf.rss_children[i] == rt and n < pf.pop_size_d):
          pf.p_temp[n,:] = pf.p_children[i,:]
          pf.rss_temp[n] = pf.rss_children[i]
          n = n + 1
          
      # UPDATE RSS AND PS
      pf.rss[1:pf.pop_size_d] = numpy.copy(pf.rss_temp[1:pf.pop_size_d])  
      pf.ps[1:pf.pop_size_d] = numpy.copy(pf.p_temp[1:pf.pop_size_d])  


      # RUN EXTINCTION EVENT
      if((gen % pf.extinction_frequency) == 0 and gen > 0 and gen < pf.generations-1):
        pf.stage = 'Extinction event'
        pf.extinction()

      # RUN ENHANCE EVENT
      if(((gen + 1) % pf.enhance_frequency) == 0):
        pf.stage = 'Enhancing Top 10% With Steepest Descent'
        pf.enhance()
        
      # SAVE POTENTIAL
      pf.output(c, gen+1, ct_id)
        
        
    # Run end function
    pf.stage = 'End of Cycle' 
    results = pf.end()   
    pf.progress = 100
    pf.time_remaining = 0.0
    display.pot_fit_show(results) 
    time.sleep(3)
    
    # Save result
    pf.result_ps = numpy.copy(pf.ps[-5,:])
    pf.result_rss = pf.rss[-5]
    pf.update_potential(pf.ps[-5,:])
    pf.cycle_ps[c_id + 1,:] = pf.ps[-5,:]
    
    potential.save_potential('fitting_' + str(c_id))
        
    

        
  def rss_threshold(f=50):
    rss_sorted = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    t = int(numpy.floor(len(rss_sorted) * (f / 100)))
    return rss_sorted[t]
    
    
  def check_improvement(params, rss): 
    # Increment Improvement Counter
    pf.since_improvement = pf.since_improvement + 1 
    

    # Update best (if better)
    if(pf.rss[-5] < 0 or rss < pf.rss[-5]):
      pf.ps[-5,:] = params
      pf.rss[-5] = rss  
      pf.improvement_time = time.time()  
      pf.since_improvement = 0
      pf.best_rss = pf.rss[-5]
    
    pf.last_rss = rss 
    
    # Save to print in display
    pf.last_ps = params
    pf.best_ps = pf.ps[-5,:]
    pf.density_list = rescale_density.max_densities()
    display.pot_fit_show()
    
    
          
  def breed(pa, pb, ca, cb):    
    state = True    
    for i in range(pf.p_width):
      rn = random.random()
      if(rn <= pf.breed_switch):
        if(state):
          state = False
        else:
          state = True
      if(state):
        pf.p_children[ca, i] = pf.ps[pa,i]
        pf.p_children[cb, i] = pf.ps[pb,i]
      else:      
        pf.p_children[cb, i] = pf.ps[pa,i]
        pf.p_children[ca, i] = pf.ps[pb,i]
        
    # Mutation (possibly)
    pf.p_children[ca, :] = pf.mutate(pf.p_children[ca, :], chance=5) 
    pf.p_children[cb, :] = pf.mutate(pf.p_children[cb, :], chance=5)    

    
    # No clones - Child A
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.p_children[ca, :] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.p_children[ca, :] = pf.p_children[ca, :] * pf.ncv * (0.5 - r)
    
    # No clones - Child B
    vary = False
    for p in range(pf.pop_size_d):
      if((pf.p_children[cb, :] == pf.ps[p,:]).all()):
        vary = True
        break
    if(vary):
      r = numpy.random.rand(pf.p_width)
      pf.p_children[cb, :] = pf.p_children[cb, :] * pf.ncv * (0.5 - r)

    
  def extinction():
    # Extinction for those below rss_thresh
    # Replace with variations of the best 25% parameters
  
    pf.extinction_counter = pf.extinction_counter + 1
    rss_thresh = pf.rss_threshold(pf.extinction_percentage)
    pf.extinction_threshold_k = rss_thresh
    
    top_size = int(pf.pop_size_d / 4)    
    pf.rss_sort = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    top_rss_min = pf.rss_sort[top_size-1]
    pf.extinction_threshold_t = top_rss_min
    
    n = 0
    for i in range(pf.pop_size_d):
      if(pf.rss[i] <= top_rss_min and n < top_size):
        pf.p_top[n,:] = pf.ps[i,:]
        n = n + 1
      if(n == top_size):
        break
    
    #pf.rss_sort = numpy.zeros((pf.pop_size_d,),)

    
    for p in range(pf.pop_size_d):
      if(pf.rss[p]>rss_thresh):   
        r_top = numpy.random.randint(top_size)
        
        pf.ps[p,:] = pf.random_p(pf.p_top[r_top], 0.01)

        pf.update_potential(pf.ps[p,:])
        pf.rss[p] = pf.get_rss() 
        pf.check_improvement(pf.ps[p,:], pf.rss[p])
        

    
  def enhance():
    # Enhance best 10% with gradient descent
    top_size = int(pf.pop_size_d / 10)    
    pf.rss_sort = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    top_rss = pf.rss_sort[top_size-1]
    
    for p in range(pf.pop_size_d):
      if(pf.rss[p] <= top_rss):
        params = gd.opt(pf.gd_rss, pf.ps[p])
        if(gd.rss_out < pf.rss[p]):
          pf.rss[p] = gd.rss_out
          pf.ps[p, :] = numpy.copy(params)
          
    
  def gd_rss(params):
    pf.update_potential(params)
    return pf.get_rss() 

    
  def estimate_calc_counter():
  
    e = 0
    
    for c in range(pf.spline_cycles + 1):
      e = e + pf.pop_size
      e = e + pf.pop_size
            
      gens = pf.generations
      if(c > 1):    
        gens = pf.spline_generations
      for gen in range(gens):
        e = e + 2 * pf.pop_size
        e = e + 2 * pf.fresh_size
        if((gen % pf.extinction_frequency) == 0 and gen > 0 and gen < pf.generations-1):
          e = e + pf.pop_size
        if((gen % pf.enhance_frequency) == 0):
          e = e + (300 * (pf.pop_size // 4))
    
    pf.calc_counter_expected = e
  
  
   
  def random_p(c=0.0, m=1.0):
    # 
    lower = pf.p_var[0,:]
    upper = pf.p_var[1,:]
    range = upper - lower
    
    # If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      c = lower + 0.5 * range 
    
    # Multiply range
    m_range = m * range
    
    # Get random parameters 0 to 1
    r = numpy.random.rand(pf.p_width)
    
    # New parameters
    p_new = c + (r - 0.5) * m_range
    
    # Return
    return p_new
    
      
    
  def mutate(params, chance=5):
    mutant = pf.random_p(params, 1.0)
    for i in range(len(params)):
      if(random.randrange(0, 100) <= chance):
        params[i] = mutant[i]
    return params
    
    
    
  def end():
    pf.update_potential(pf.ps[-5,:])
    rss = pf.get_rss() 
    
    results = []
    results.append("RSS: " + str(rss))
    
  
    for bp_id in range(bp.bp_configs_count):
      results.append("alat: " + str(bp.calc_alat[bp_id]))
      results.append("v0: " + str(bp.calc_v0[bp_id]))
      results.append("e0: " + str(bp.calc_e0[bp_id]))
      results.append("b0: " + str(bp.calc_b0[bp_id]))
      results.append("b0/GPA: " + str(160.230732254e0 * bp.calc_b0[bp_id]))
      return results
  
  
  # To Do
  def sort_parameters():
    pf.rss[0:pf.pop_size_d] = sort.sort_1d_dp_asc(pf.rss[0:pf.pop_size_d])
    
    #kt = sort.key_table
    #pf.ps[0:pf.pop_size_d] = sort.match_key_table_2d_dp_kt(pf.ps[0:pf.pop_size_d], kt)
    
  
  def output(cycle, gen, ct_id):
  
    pot_dir = g.dirs['fitting'] + "/fits"
    std.make_dir(pot_dir)   
  
    pf.output_counter = pf.output_counter + 1
    pf.update_potential(pf.ps[-5,:])
  
    file_counter = str(pf.output_counter)    
    while(len(file_counter)<4):
      file_counter = "0" + file_counter
      
    cycle_counter = str(cycle)    
    while(len(cycle_counter)<4):
      cycle_counter = "0" + cycle_counter  
      
    gen_counter = str(gen)    
    while(len(gen_counter)<4):
      gen_counter = "0" + gen_counter 
      
      
    # Make dir
    fit_dir = pot_dir + "/" + file_counter + "_" + str(ct_id) + "_" + cycle_counter + "_" + gen_counter  
    std.make_dir(fit_dir)       
    
    # Output plots
    potential.plot_python_potentials(fit_dir)
    potential.plot_fortran_potentials(fit_dir)
    
    # Output file
    fh = open(fit_dir + "/rss.txt", 'w')
    fh.write("RSS: " + str(pf.rss[-5]) + "\n")   
    fh.close()
    
    
    fh = open(fit_dir + "/parameters.txt", 'w')
    fh.write("Cycle:  " + cycle_counter + "\n")   
    fh.write("Gen:    " + gen_counter + "\n")   
    fh.write("Parameters" + "\n") 
    for i in range(len(pf.ps[-5,:])):
      n = str(i)
      while(len(n)<3):
        n = '0' + n
      fh.write(n + "   " + str(pf.ps[-5,i]) + "\n")  
    fh.write("RSS: " + str(pf.rss[-5]) + "\n")    
    fh.close()
        
    # Write log of all end of gen best rss    
    fh = open(g.dirs['fitting'] + "/rss.txt", 'a')
    fh.write("Cycle:  " + cycle_counter + "\n")   
    fh.write("Gen:    " + gen_counter + "\n")   
    fh.write("RSS:    " + str(pf.rss[-5]) + "\n")  
    fh.write("" + "\n")    
    fh.close()
    
    
    potential.save_potential(fit_dir + "/pot", "potential")
    
    
    
    