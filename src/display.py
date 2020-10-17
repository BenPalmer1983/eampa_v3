import numpy
import os
from globals import g
from pot_fit import pf

class display:

  last_stage = None

  def clear():    
    os.system('cls' if os.name == 'nt' else 'clear') 
  
  def print_line(w=140):
    for i in range(w):
      print("#", end="")
    print()

  def output():
    try:
      output = int(g.inp['display']['output'])
    except:
      output = 1    
    if(output == 1):
      display.output_1()
    if(output == 2):
      display.output_2()
    if(output == 3):
      display.output_3()
    
    
    
  def output_1():
    if(display.last_stage == None or display.last_stage != g.pfdata['stage']):
      display.last_stage = g.pfdata['stage']
      print("Stage: " + g.pfdata['stage'])
  
    print("   " , g.pfdata['rss']['counter'], g.pfdata['rss']['current'], g.pfdata['rss']['best'])
  
    
    
  def output_2():
    if(display.last_stage == None or display.last_stage != g.pfdata['stage']):
      display.last_stage = g.pfdata['stage']
      print("Stage: " + g.pfdata['stage'])
  
    print("   " , g.pfdata['rss']['counter'], g.pfdata['rss']['current'], g.pfdata['rss']['best'])
  
  
  def output_3(results=None):

    # HEADER
    display.header_3()

    # PRINT RSS
    print("# RSS:                ", g.pfdata['rss']['current'], "  [",g.pfdata['rss']['best'],"]") 
    display.print_line()

    # Best BP
    bp_best = g.pfdata['bp_best']
    
    # PRINT VALUES
    print("#        ID   a0       e0       B0       C11      C12      C44 ")
    for bp_id in g.bp_results['input'].keys():
      print(  display.pad_r("# Calc: ", 9)
            + display.pad_r(bp_id,4) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['a0'],8) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['e0'],8) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['b0_gpa'],8) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
            + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    for bp_id in bp_best['input'].keys():
      print(  display.pad_r("# Best: ", 9)
            + display.pad_r(bp_id,4) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['a0'],8) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['e0'],8) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['b0_gpa'],8) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
            + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    for bp_id in g.bp_known['input'].keys():      
      print(  display.pad_r("# Known: ", 9)
            + display.pad_r(bp_id,4) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['a0'],8) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['e0'],8) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['b0_gpa'],8) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
            + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    display.print_line()
    
    
    # PRINT PARAMETERS
    
    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("Fn: " + str(fn) + "   Type: spline")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['fit_parameters'][0,i],8), end='')
          print()
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        print("Fn:" + str(fn) + "[A] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['a_params'][i],8), end='')
        print()
    display.print_line()
    
    
    
    
  def finish():
    try:
      option = int(g.inp['display']['output'])
    except:
      option = 1    
    if(option == 1):
      display.finish_1()
    if(option == 2):
      display.finish_2()
    if(option == 3):
      display.finish_3()
    
    
  def finish_3(results=None):
    # HEADER
    display.header_3()
    
    # PRINT RSS
    display.print_line()
    print("# BEST RSS:                  ",g.pfdata['rss']['best'],"") 
    display.print_line()
  
    
    
    
    """

    print("###################################################################")
    print("#  Pop/Fresh size  :   " + display.pad_l(pf.pop_size)  + " /  " +  display.pad_l(pf.fresh_size))
    print("#  Generations:        " + display.pad_l(pf.generations) )
    print("#  Spline Cycles:      " + display.pad_l(pf.spline_cycles) + " /  " + display.pad_l(pf.spline_generations) + "")
    print("#  Ext/Enhance Freq:   " + display.pad_l(pf.extinction_frequency) + " /  " + display.pad_l(pf.enhance_frequency) + "")    
    print("###################################################################")
    print("# Timer:              " + display.pad_l(time.time() - pf.start_time, 10) + " (Estimate: " + display.pad_l(pf.time_estimate, 10) + ")")
    print("# Stage:              " + pf.stage) 
    print("# Cycle:              " + display.pad_l(pf.this_cycle, 10))  
    print("###################################################################")
    print("# Gen:                " + display.pad_l(pf.this_gen, 10))  
    print("# Exctinctions:       " + display.pad_l(pf.extinction_counter, 10) 
          + " (" + display.pad_l(pf.extinction_threshold_t, 10) + "/" 
          + "" + display.pad_l(pf.extinction_threshold_k, 10) + ")") 
    print("# Since Improvement:  " + display.pad_l(pf.since_improvement, 10))  
    print("#    ")  
    print("# RSS:           ", pf.last_rss, "  [",pf.best_rss,"]") 
    print("#    ")   
    for bp_id in range(len(pf.bulk_properties)):
      print("#                     a0  " + display.pad_l(pf.bulk_properties[bp_id,3], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,6], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,0], 10) + " ") 
      print("#                     e0  " + display.pad_l(pf.bulk_properties[bp_id,4], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,7], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,1], 10) + " ") 
      print("#                     b0  " + display.pad_l(pf.bulk_properties[bp_id,5], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,8], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,2], 10) + " ") 
      ec_list = [11,22,33,44,55,66,12,13,23]
      for ec in range(len(ec_list)):
        print("#                     c" + str(ec_list[ec]) + " " + display.pad_l(pf.bulk_properties[bp_id,ec + 21], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,ec + 31], 10) + "   " + display.pad_l(pf.bulk_properties[bp_id,ec+11], 11) + " ") 
      
     
     
     
     
    print("#    ")  
    print("# display:   ", display.bar(pf.display, 25))  
    print("# TIME LEFT:  ", round(pf.time_remaining,2))  
    print("#   ")  
    print("# Parameters (this): ", end="")
    if(type(pf.last_ps) == numpy.ndarray):
      for pn in range(len(pf.last_ps)):
        print(pf.last_ps[pn], end=" ")
        if((pn+1) % 5 == 0):
          print()
          print("#                   ", end="")
    print("")
    print("# Parameters (best): ", end="")
    if(type(pf.best_ps) == numpy.ndarray):
      for pn in range(len(pf.best_ps)):
        print(pf.best_ps[pn], end=" ")
        if((pn+1) % 5 == 0):
          print()
          print("#                   ", end="")
    print("")
    print("#   ")  
    print("# Density List: ", pf.density_list)  
    print("#   ")  
    
    
    #    pf.last_ps = params
    #pf.best_ps = params
    
    print("###################################################################")
    
    if(type(results) == list):
      for r in results:
        print("# " + str(r))  
      print("###################################################################")
        
    
    #for bp_id in range(bp.bp_configs_count):
    
    #pf.since_improvement

    """
 
  def header_3():
    # CLEAR
    display.clear()  
  
  
    # Extinction
    e = (g.fit['exct_every'] - g.pfdata['extinction']['counter']) - 1
    if( e == "0"):
      e_print = "This generation"
    else:
      e_print = "In " + str(e) + " generations"
        
    
    # PRINT HEADER
    display.print_line()
    print("# Cycle: " + str(g.pfdata['cycle']['counter']) + " of " 
          + str(g.pfdata['cycle']['total_cycles']) + "    "
          + "Generation: " + str(g.pfdata['generation']['counter']) + " of " 
          + str(g.pfdata['generation']['total_generations']) + "    ")
    print("# Pop Size:    " + str(g.fit['pop_size']) + "      Fresh Size:    " +  str(g.fit['pop_size']))
    display.print_line()
    
    line = ['','','','','']
    # Col 1
    line[0] = display.pad_r_always("# Timer:              " + display.pad_l(time.time() - g.times['start'], 16), 60)
    line[1] = display.pad_r_always("# Stage:              " + g.pfdata['stage'], 60)
    line[2] = display.pad_r_always("# RSS Counter:        " + display.pad_l(g.pfdata['rss']['counter'], 16), 60)
    line[3] = display.pad_r_always("# Since Improvement:  " + display.pad_l(g.pfdata['rss']['since_improvement'], 16), 60)
    line[4] = display.pad_r_always("# Next Extinction:    " + e_print, 60)
    
    # Col 2
    line[0] = line[0] + "# " + display.pad_r_always("TOP 10", 21)
    
    for n in range(10):
      ln = (n%4) + 1
      if(g.pfdata['top']['filled']):
        line[ln] = line[ln] + display.pad_r_always(g.pfdata['top']['rss'][n], 20)
      else:
        tn = n + (g.pfdata['top']['size'] - g.pfdata['top']['counter']) 
        if(tn < g.pfdata['top']['size']):
          line[ln] = line[ln] + display.pad_r_always(g.pfdata['top']['rss'][tn], 20)
    
    for n in range(5):
      print(line[n])
  
    display.print_line()
      
      
      
 
 
  @staticmethod
  def pad_r(inp, p=7):
    if(inp == None):
      return ""    
    if(type(inp) == "float64"):
      inp = numpy.round(inp, p-3)
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]
    
  @staticmethod
  def pad_r_always(inp, p=7):
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]  
    
  @staticmethod
  def pad_l(inp, p=7):
    if(inp == None):
      return ""      
    out = str(inp).strip()  
    while(len(out)<p):
      out = " " + out     
    return out[0:p]
    
  @staticmethod
  def bar(p, w=25):
    out = ''
    p_num = p
    p = p * (25 / 100)
    for i in range(w):
      if(i <= p):
        out = out + '#'
      else:
        out = out + '_'
    out = out + '   ' + str(p_num) + '%'
    return out
    
 
  """
  start_time = 0
  since_improvement = 0
  this_gen = 0
  this_cycle = 0
  total_generations = 0
  last_rss = 0.0
  best_rss = 0.0
  """

