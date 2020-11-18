##########
# GLOBALS

import numpy

class g: 
  
  dirs = {
         'wd': 'wd',
         }
  
  sub_dirs = { 
         'input': 'input', 
         'output': 'output',   
         'results': 'results',   
         'plots': 'plots',     
         'configs': 'configs',  
         'fitting': 'fitting', 
         'fitting_generations': 'fitting/generations', 
         'fitting_finished': 'fitting/finished', 
         }
  
  times = {
          'start' : 0.0,
          'end' : 0.0,
          'duration' : 0.0,
          }
          
  pot_functions = {
                  'pot_name': '',
                  'pot_dir': '',
                  'zbl_file': '',
                  'functions': [],
                  'zbl': [],
                  'functions_original': [],
                  }
                  
  configs = {
            'config_files': [],
            'configs': [],
            }  
            
  mask = {}
            
      

  dft_energy_adjustments = {}  
  bulk_properties = []
  bp_ids = {}  
  labels = {}
  grouplabels = {}
  groupelement = {}
  
  # Read in from input
  run_type = 'efs'
  wd_type = {}
  rss_weights = {}
  fit = {}
  fit_results = {}


  ######################################
  # Calculated results
  # rss, fit
  ######################################
  
  efs_known = {'ok': False,}
  efs_results = {'ok': False,}
  efs_results_best = {'ok': False,}
  
  bp_known = {'ok': False,}
  bp_results = {'ok': False,}
  bp_results_best = {'ok': False,}
  
  benchmark = {'configs': 0, 'total_atoms': 0, 'total_interactions': 0, 'total_time': 0.0, 'configspersec': 0.0, 'atomspersec': 0.0, 'interationspersec': 0.0,}
  

  #RSS
  rss = {'current': None, 'best': None, 'counter': None, 'since_improvement': None, 'log': [], 'efs': {'ok': False, 'cc': 0,}, 'bp': {'ok': False, 'cc': 0,}}
  
  ######################################
  # Potential Fitting
  ######################################
  
  pfdata = {}
  
  


  ######################################
  # Other settings
  ######################################


  tab_size = 1001
  tab_width = 4
  outputs = True
  results_fh = None
  log_fh = None         
  file_counter = 0 
         
  def file_name():
    globals.file_counter = globals.file_counter + 1
    name = "file_"
    file_counter_str = str(globals.file_counter)
    while(len(file_counter_str) < 6):
      file_counter_str = '0' + file_counter_str
    name = name + file_counter_str    
    return name
         
         