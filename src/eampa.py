######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
#from tendl import tendl
#from isotopes import isotopes
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from read_input import read_input
from memory import memory
from labels import labels
from potential import potential
from potential_vary import potential_vary
from configs import configs
from b_props import b_props
from efs_calc import efs_calc
from bp_calc import bp_calc
from es_calc import es_calc
from rss_calc import rss_calc
from pot_fit import pf
from trial import trial
from relax_calc import relax_calc



from eampa_lib.f_es import es
from eampa_lib.f_efs import efs
from eampa_lib.f_bp import bp
from eampa_lib.f_sorting import sort
from eampa_lib.f_interp import interp
from eampa_lib.f_spline import spline
from eampa_lib.f_fnc import fnc
from eampa_lib.f_bp import polyfit
from eampa_lib.f_relax import relax



class eampa:
 
  def run():
    print("RUNNING")
    
    # Copy Input
    eampa.copy_input()   
    
    # Read Input
    read_input.run()    
    
    # Set memory
    memory.run() 
        
    # Load potentials
    potential.load()
    
    # Load configs
    configs.load()
        
    # Bulk Properties
    b_props.load()
    
    # Surface Energy etc
    es_calc.load()
    
    # Convert to ev/ang etc and adjust energies
    configs.complete()

    labels.output()    
    configs.output()
    
    print("RUN TYPE: " + g.run_type)
    time.sleep(0.4)
    
    
    if(g.run_type == 'e'):
      efs_calc.run_energy()
    elif(g.run_type == 'ef'):
      efs_calc.run_energy_force()
    elif(g.run_type == 'efs'):
      efs_calc.run_energy_force_stress()
    elif(g.run_type == 'bp'):
      bp_calc.run()
    elif(g.run_type == 'es'):
      es_calc.run()
    elif(g.run_type == 'rss'):
      rss_calc.run()
    elif(g.run_type == 'fit'): 
      pf.run()
    elif(g.run_type == 'plot'): 
      potential.run()
    elif(g.run_type == 'trial'): 
      trial.run()
    elif(g.run_type == 'relax'): 
      relax_calc.run()
    else: 
      trial.run()
      
   
  def copy_input():
    try: 
      std.copy(g.inp['potential']['dir'], g.dirs['input'])
    except:
      pass
    try: 
      std.copy(g.inp['configs']['dir'], g.dirs['input'])
    except:
      pass
    try: 
      std.copy(g.inp['dft']['dir'], g.dirs['input'])
    except:
      pass
    try: 
      std.copy(g.inp['bp']['dir'], g.dirs['input'])
    except:
      pass
    
    
  """  
  def run_type():
    # Types
    # config       just calculate energy/forces/stress of configs
    # bp           calculate bulk properties (bulk modulus, elastic constants etc)
    # rss          
    # fit  
  
    # Default
    g.run_type = g.inp['run']['type'].lower().strip()
  """
  
  
  

  def exit():
    print("End of Program")
    exit()
  