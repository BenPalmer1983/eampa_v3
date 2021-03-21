################################################################
#    Main Program
#
#
#
#
################################################################


#!/bin/python3
########################################################################
import os
import time
import datetime
import re
import sys
import shutil
import numpy
import matplotlib.pyplot as plt
import copy
import random
import hashlib
from eampa_lib.f_es import es
from eampa_lib.f_efs import efs
from eampa_lib.f_bp import bp
from eampa_lib.f_sorting import sort
from eampa_lib.f_interp import interp
from eampa_lib.f_spline import spline
from eampa_lib.f_fnc import fnc
from eampa_lib.f_bp import polyfit
from eampa_lib.f_relax import relax

###########################################
#  CLASS
###########################################
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
  pot_labels = []

  pot_fn = {}
                  
  configs = {
            'config_files': [],
            'configs': [],
            'config_results': [],
            }  
            
  mask = {}
            
  dft_energy_adjustments = {}  
  bulk_properties = []
  bp_ids = {}  
  labels = {}
  fgroups = {}
  grouplabels = {}
  groupelement = {}
  
# Read in from input
  run_type = 'efs'
  wd_type = {}
  rss_weights = {}
  rss_max_density = {}
  fit = {}
  fit_results = {}
  
# Top 100 parameters/rss
  top_parameters = []

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
  rss = {'current': None, 'best': None, 'counter': None, 'since_improvement': None, 'log': [], 'efs': {'ok': False, 'cc': 0,}, 'bp': {'ok': False, 'cc': 0,}, 'residual': []}
  
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
         
###########################################
#  CLASS st
###########################################
class std:

  @staticmethod
  def file_to_list(file_name, clean=False):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(clean):
        line = line.strip()
        if(line != ""):
          file_data.append(line)          
      else:
        file_data.append(line[0:-1])
# Return
    return file_data
    
  @staticmethod
  def split_fields(line, sep=" "):
    out = line.split(sep)
    key = out[0]
    value = out[1]
    value_out = ''    
    indata = False
    for char in value:
      if(indata and char != '"'):
        value_out = value_out + char
      elif(indata and char == '"'):
        indata = False
      elif(not indata and char == '"'):
        indata = True
    return key, value_out
    
  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out   
    
  @staticmethod
  def to_fields(line, sep=" "):
    out = []
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == sep and last_char == sep)):
        if(char == sep):
          temp = temp.strip()
          if(temp != ""):
            out.append(temp)
            temp = ''
        else:
          temp = temp + char
    
    temp = temp.strip()
    if(temp != ""):
      out.append(temp)      
    return out    
    
  @staticmethod
  def make_dir(dir):
    dirs = dir.split("/")
    try:
      dir = ''
      for i in range(len(dirs)):
        dir = dir + dirs[i]
        if(not os.path.exists(dir) and dir.strip() != ''):
          os.mkdir(dir) 
        dir = dir + '/'
      return True
    except:
      return False
      
  @staticmethod
  def copy(src, dest):  
    try:
      std.make_dir(dest)
      src_f = src.split("/")
      if(os.path.isdir(src)):
        shutil.copytree(src, dest + '/' + src_f[-1])      
      else:
        shutil.copyfile(src, dest + '/' + src_f[-1])
      return True
    except:
      return False
      
  @staticmethod
  def path(dir, file):  
    dir = dir.strip()
    file = file.strip()
    if(dir == ''):
      return file
    else:
      if(dir[-1] == '/'):
        return dir + file
      else:      
        return dir + '/' + file
    
  @staticmethod
  def remove_comments(content):
    data = ''
    i = 0
    for line in content:
      if(i > 0):
        data += '\n'
      data += line
      i = i + 1
    out = ''
    indata = 0
    incomment = 0
    for i in range(len(data)):
# Get char and next char
      char = data[i]
      next = None
      prev = None
      if(i < len(data)-1):
        next = data[i + 1]
      if(i > 0):
        prev = data[i - 1]
# If in '  '
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
# If in "  "
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0):
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 1 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "!"):
          incomment = 2
        elif(incomment == 2 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
    return out.split("\n")    
    
# Remove comments from a block of data/text
  @staticmethod
  def remove_comments_data(data):
    out = ""
    n = 0
    inquotes = 0
    incomment = 0
    while n < len(data):
# Get char and next char
      char = data[n]
      next = None
      prev = None
      if(n < len(data)-1):
        next = data[n + 1]
      if(n > 0):
        prev = data[n - 1]
        
# If in '  '
      if(inquotes == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(inquotes == 1 and char == "'" and last_char != "\\"):
        out = out + char
        inquotes = 0
# If in "  "
      elif(inquotes == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(inquotes == 2 and char == '"' and last_char != "\\"):
        out = out + char
        inquotes = 0
# If not inside quotes
      elif(inquotes == 0):
# Comment on a line
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 0 and char == "!"):
          incomment = 1
        elif(incomment == 0 and char == "#"):
          incomment = 1    
# Comment on line close
        elif(incomment == 1 and char == "\n"):
          incomment = 0
# Comment block
        elif(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
# Increment counter
      n = n + 1
    return out        

# Single spaces, tabs to spaces
  @staticmethod
  def prep_data(content):
    out = []
    for line in content:
      line_new = std.prep_data_line(line)
      if(line_new != ''):
        out.append(line_new)
    return out  
      
  @staticmethod
  def prep_data_line(line): 
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(char == '\t'):
        char = ' '
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == ' ' and last_char == ' ')):
        temp = temp + char       
      last_char = char  
    return temp.strip()    
    
  @staticmethod
  def remove_quotes(inp): 
    if(isinstance(inp, list)):    
      for i in range(len(inp)):
        inp[i] = std.remove_quotes(inp[i])        
      return inp
    else:
      inp = inp.strip()
      if(inp[0] == '"' and inp[-1] == '"'):
        return inp[1:-1]
      if(inp[0] == "'" and inp[-1] == "'"):
        return inp[1:-1]
      return inp
      
  @staticmethod
  def config_file_to_list(file_name):
# Init variable
    file_data = []
# Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(line.strip() != ""):
        line = line.strip()
        line = std.remove_comments(line)
        line = std.prep_data_line(line)
        fields = std.to_fields(line)
        file_data.append(fields)         
# Return
    file_data = std.remove_quotes(file_data)
    return file_data
    
  @staticmethod
  def get_dir(file_path):
    directory = ''
    read = False
    for i in range(len(file_path)):
      if(read):
        directory = file_path[-1-i] + directory
      if(file_path[-1-i] == "/"):
        read = True
    return directory
  
  @staticmethod
  def read_csv(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
# Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line
      fh.close()
# Remove comments
      file_data = std.remove_comments_data(file_data)
# Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
    return data
     
  @staticmethod
  def read_csv_array(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
# Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line.strip() + '\n'
      fh.close()
# Remove double spaces
      file_data = std.one_space(file_data)
# Remove comments
      file_data = std.remove_comments_data(file_data)
# Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
          
      lst = []
      for i in range(len(data)):
        row = []
        for j in range(len(data[0])):
          try:
            row.append(float(data[i][j]))
          except:
            pass
        lst.append(row)
        
      arr = numpy.zeros((len(lst),len(lst[0]),),)
      for i in range(len(lst)):
        for j in range(len(lst[0])):
          arr[i,j] = float(lst[i][j])
      return arr
    return None
    
  @staticmethod
  def write_csv(filename, arr):  
    fh = open(filename, 'w')
    for i in range(len(arr)):
      for j in range(len(arr[i])):
        fh.write(std.float_padded(arr[i,j],8) + " ")
      fh.write("\n")
    fh.close()
  
  @staticmethod
  def option(input):
    input = input.strip().upper()
    if(input[0:1] == "Y"):
      return True
    elif(input[0:2] == "ON"):
      return True
    elif(input[0:1] == "T"):
      return True
    else:
      return False
    
  @staticmethod
  def float_padded(inp, pad=7):
    out = float(inp)
    out = round(out, pad-3)
    out = str(out)  
    while(len(out)<pad):
      out = out + " "      
    return out[0:pad]
    
  @staticmethod
  def write_file_line(fh, title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    fh.write(line)  
  
  @staticmethod
  def print_file_line(title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    print(line,end='')  
    
  def pad(inp, width):
    out = str(inp)
    while(len(out)<width):
      out = out + " "
    return out
    
  def mem_value(strin):
    num = '0123456789.'
    strin = strin.upper().strip()
    value = ''
    unit = ''
    for c in strin:
      if(c not in num):
        break
      else:
        value = value + c
    for c in strin:
      if(c not in num):
        unit = unit + c
     
    value = float(value)
    if(unit == 'B'):
      return value
    if(unit == 'KB'):
      return 1000 * value
    if(unit == 'MB'):
      return 1000000 * value
    if(unit == 'GB'):
      return 1000000000 * value
    
  @staticmethod
  def dict_to_str(d, pre=''):
    out = ''
    for k in sorted(d):
      if(type(k) == dict or type(k) == list):
        out = out + '\n' + std.dict_to_str(k, pre + '  ')      
      else:
        out = out + pre + str(k) + ': ' + str(d[k]) + '\n'
    return out
    
###########################################
#  CLASS read_confi
###########################################
class read_config:
  
  @staticmethod
  def read_file(file_path):
  
# Input dictionary
    input = {}
  
# READ DATA
    d = []
    fh = open(file_path, 'r')
    for line in fh:
      line = line.strip()
      if(len(line) > 0 and line[0] != "#"):
        d.append(line)      
    fh.close()
    
# Count commands
    commands = {}
    for line in d:
      fields = read_config.split_by(line, ' ')
      c = fields[0].lower()
      if(c in commands.keys()):
        commands[c] = commands[c] + 1
      else:
        commands[c] = 1
        
# Prepare input dictionary
    for k in commands.keys():
      if(commands[k] == 1):
        input[k] = None
      else:
        input[k] = []
    
# Read Data into input
    for line in d:
      fields = read_config.split_by(line, ' ')
      fkey = fields[0].lower()
      
      fd_size = {}
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")
        fc = fs[0].lower()
        if(fc in fd_size.keys()):
          fd_size[fc] = fd_size[fc] + 1
        else:
          fd_size[fc] = 1
          
# Prepare dictionary
      fd = {} 
      for k in fd_size.keys():
        if(fd_size[k] == 1):
          fd[k] = None
        else:
          fd[k] = []        
        
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")     
        fc = fs[0].lower()        
        fs = read_config.split_by(fs[1], ',')         
        fs = read_config.store(fs)
        
        if(fd_size[fc] == 1):
          if(len(fs) == 1):
            fd[fc] = read_config.store(fs[0])
          else:
            fd[fc] = read_config.store(fs)
        else:
          if(len(fs) == 1):
            fd[fc].append(read_config.store(fs[0]))
          else:
            fd[fc].append(read_config.store(fs))
            
      if(commands[fkey] == 1):
        input[fkey] = fd
      else:
        input[fkey].append(fd)  

    return input
        
  @staticmethod  
  def split_by(line, sep=' ', ignore_double_sep=True):
    last_char = None
    in_quotes = 0
    fields = []
    temp_line = ""
    
    for char in line:
      if(char == "'" and in_quotes == 0 and last_char != "\\"):
        in_quotes = 1
      elif(char == "'" and in_quotes == 1 and last_char != "\\"):
        in_quotes = 0
      elif(char == '"' and in_quotes == 0 and last_char != "\\"):
        in_quotes = 2
      elif(char == '"' and in_quotes == 2 and last_char != "\\"):
        in_quotes = 0
      elif(in_quotes > 0):
        temp_line = temp_line + char
      elif(in_quotes == 0 and char != sep):
        temp_line = temp_line + char
      elif(char == sep and last_char == sep and ignore_double_sep):
        pass
      elif(char == sep):
        fields.append(temp_line)
        temp_line = "" 
    if(temp_line != ""):
      fields.append(temp_line)
    
    return fields
    
  @staticmethod
  def store(inp):  
    if(isinstance(inp, list)):
      for i in range(len(inp)):
        try:
          if('.' in inp[i]  or 'e' in inp[i]):
            inp[i] = float(inp[i])
          else:
            inp[i] = int(inp[i])
        except:
          pass
    else:
      try:
        if('.' in inp or 'e' in inp):
          inp = float(inp)
        else:
          inp = int(inp)
      except:
        pass
    return inp
      
###########################################
#  CLASS eamp
###########################################
class eampa:
 
  def run():
    print("RUNNING")
        
# Read Input
    read_input.run()    
    
# Setup Dirs
    setup_dirs.run()
    
# Copy Input
#eampa.copy_input()
    
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
  
###########################################
#  CLASS read_inpu
###########################################
class read_input:
 
  def run():
    main.log_title("Read Input")
  
    read_input.run_type()
    read_input.wd()
    read_input.rss_weights()
    read_input.rss_max_density()
    read_input.fitting()
    read_input.fit()
    read_input.fit_results()
    read_input.bp()
    read_input.mask()
    read_input.dft_ea()

# READ TYPE
  def run_type():
  
# DEFAULT
    g.run_type = 'efs'
    
# TRY READING
    try:
      g.run_type = g.inp['run']['type'].lower().strip()
    except:
      pass
      
# SAVE
    main.log(g.run_type)

# READ TYPE
  def wd():
    g.wd_type = {
    'option': 1,
    }    
# TRY READING
    for k in g.wd_type.keys():
      try:
        g.wd_type[k] = float(g.inp['wd'][k])
      except:
        pass

# READ
  def rss_weights():
  
# DEFAULT
    g.rss_weights = {
    'config': 1.0,
    'force': 1.0,
    'energy': 1.0,
    'stress': 1.0,
    'bp': 1.0,
    'a0': 1.0,
    'e0': 1.0,
    'b0': 1.0,
    'ec': 1.0,
    'g': 1.0,
    'e': 1.0,
    'v': 1.0,
    'negec': 1.0,
    }
    
# TRY READING
    for k in g.rss_weights.keys():
      try:
        g.rss_weights[k] = float(g.inp['rss_weights'][k])
      except:
        pass
# READ
      
# SAVE
    main.log(std.dict_to_str(g.rss_weights))
      
# READ
  def rss_max_density():   
      
# DEFAULT
    g.rss_max_density = {          
    'min': 0.2,
    'max': 0.8,
    'scale_factor': 10.0,
    'scale_exponent': 4.0,
    'zero_density_factor': 1.0e8,
    }
    
# TRY READING
    for k in g.rss_max_density.keys():
      try:
        g.rss_max_density[k] = float(g.inp['rss_max_density'][k])
      except:
        pass
# READ
      
# SAVE
    main.log(std.dict_to_str(g.rss_max_density))
    
# READ
  def fitting():   
      
# DEFAULT
    g.fitting = {
    'oversized_parameters': [10.0,0.05,0.05,0.05],   
    'top_parameters': 100,
    'load_top_parameters': 10,
    }
    if('fitting' in g.inp.keys()):
      for k in g.fitting.keys():
        try:
          g.fitting[k] = g.inp['fitting'][k]
        except:
          pass

  def fit():
  
# List of fit types
    g.fit = []
    if('fit' in g.inp.keys()):
      g.fit.append(read_input.read_fit('fit'))
    for i in range(100):
      if("fit" + str(i) in g.inp.keys()):
        g.fit.append(read_input.read_fit("fit" + str(i)))

  def read_fit(inp_key):  
    
    fit_data = {
    'type': None,
    'random_size': 100,
    'cycles': 0,
    'gens': 0,
    'spline_cycles': 0,
    'spline_gens': 0,
    'pop_size': 20,
    'fresh_size': 10,
    'gen_variation_multiplier': 0.1,
    'exct_factor': 0.5,
    'exct_every': 5,
    'exct_var': 0.1,
    'exct_top_bias': 0.5,
    'rescale_density': 0,
    'wide_start': 0.5,
    'wide_end': 10.0,
    'mutate_chance': 0.01,
    'mutate_scale': 1.0,
    'no_clones': True,
    'no_clone_var': 0.05,
    'fresh_ws': 0.1,
    'fresh_we': 1.5,
    'enhance_every': 10,
    'enhance_top': 5,
    'gen_var_factor': 1.0,
    'pool_size': 1000,
    'sane_seeds_a': 50,
    'sane_seeds_b': 200,
    'sane_fraction': 0.5, 
    'sa_loops_t': 10,
    'sa_loops_i': 100,
    'sa_temp_start': 10,
    'sa_temp_end': 0.1,
    'sa_step': 0.01,
    'sa_step_factor': 0.3,
    }

# TRY READING
    for k in fit_data.keys():
      try:
        fit_data[k] = g.inp[inp_key][k]
      except:
        pass
        
# POP SIZE - must be even
    if(fit_data['pop_size'] < 2):
      fit_data['pop_size'] = 2
    if(fit_data['pop_size'] % 2 != 0):
      fit_data['pop_size'] = fit_data['pop_size'] + 1
      
# FRESH SIZE - must be even
    if(fit_data['fresh_size'] < 2):
      fit_data['fresh_size'] = 2
    if(fit_data['fresh_size'] % 2 != 0):
      fit_data['fresh_size'] = fit_data['fresh_size'] + 1

    if(fit_data['type'] is None):
      fit_data['type'] = 'RANDOM'

    return fit_data

  def fit_results():
      
# DEFAULT
    g.fit_results = {
    'results_dir': 'results',
    'pot_dir': None,
    'pot_name': None,
    }
    
# TRY READING
    for k in g.fit_results.keys():
      try:
        g.fit_results[k] = g.inp['fit_results'][k]
      except:
        pass
      
# SAVE
    main.log(std.dict_to_str(g.fit_results))
    
  def bp():
      
# DEFAULT
    g.bp_input = {
    'dir': '',
    'bp_file': None,
    'eos_size': 10,
    'eos_strain': 0.005,
    'ec_size': 10,
    'ec_strain': 0.005,
    }

# TRY READING
    for k in g.bp_input.keys():
      try:
        g.bp_input[k] = g.inp['bp'][k]
      except:
        pass
        
# SAVE
    main.log(std.dict_to_str(g.bp_input))

  def mask():
    g.mask = {}
    if('mask' in g.inp):
      try:
        for k in g.inp['mask']:
          g.mask[k.upper()] = g.inp['mask'][k].upper()
      except:
        pass
        
# SAVE
    main.log(std.dict_to_str(g.mask))

  def dft_ea():
  
    dftea = {}
    if('dft' in g.inp):
      try:
        for k in g.inp['dft']:
          dftea[k.upper()] = g.inp['dft'][k]
      except:
        pass
        
    for k in dftea.keys():
      label_str, label_id = labels.add(k)  
      
      atom_count = int(dftea[k][0])
      relaxed_energy = float(dftea[k][1])
      relaxed_energy_unit = str(dftea[k][2])
      coh_energy = float(dftea[k][3])
      coh_energy_unit = str(dftea[k][4])
      
      relaxed_dft_ev = units.convert(relaxed_energy_unit, "EV", relaxed_energy / atom_count) # relaxed per atom in eV
      coh_ev = units.convert(coh_energy_unit, "EV", coh_energy)
      apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
      g.dft_energy_adjustments[label_id] = {
                                            'label_id': label_id,
                                            'label_text': label_str,
                                            'atom_count': atom_count,
                                            'relaxed_energy': relaxed_energy,
                                            'relaxed_energy_unit': relaxed_energy_unit,
                                            'coh_energy': coh_energy,
                                            'coh_energy_unit': coh_energy_unit,
                                            'calc_relaxed_dft_ev': relaxed_dft_ev,
                                            'calc_coh_ev': coh_ev,
                                            'calc_apaev': apaev,
                                           }
      
# Save
    main.log(std.dict_to_str(g.dft_energy_adjustments))
  
###########################################
#  CLASS setup_dir
###########################################
class setup_dirs:
 
  def run():
    main.log_title("Setup Directories")
 
    if(g.wd_type == 1):
      
# SET UP DIRS
      now = datetime.datetime.now()
      date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
      wd_prefix = now.strftime('wd/%Y%m%d_%H%M%S')
  
# Set wd
      g.dirs['wd'] = wd_prefix
    
    else:
  
# Set wd
      g.dirs['wd'] = 'wd/wd'
      
    main.log("WD: " + g.dirs['wd'])
    
# Set sub dirs
    if('input' not in g.sub_dirs.keys()):
      g.sub_dirs['input'] = 'input'  
    for sd in g.sub_dirs.keys():
      g.dirs[sd] = g.dirs['wd'] + '/' + g.sub_dirs[sd]
  
# MAKE DIRS
    for d in g.dirs.keys():
      std.make_dir(g.dirs[d])
      main.log(str(d) + "   " + str(g.dirs[d]))
    
###########################################
#  CLASS memor
###########################################
class memory:

  def run():
  
    main.log_title("Memory")
  
    g.memory = {}

# Defaults
    g.memory['bp'] = {}    
    g.memory['efs'] = {}
      
# For BP module
    try:
      mem = str(g.inp['mem']['bp'])
    except:
      mem = "500MB"
    mem = std.mem_value(mem)    
    g.memory['bp']['mem'] = mem   
    g.memory['bp']['c'] = int(1 * (mem / 7840))
    g.memory['bp']['g'] = int(12 * (mem / 7840))
    g.memory['bp']['nl'] = int(100 * (mem / 7840))
      
# For EFS module
    try:
      mem = str(g.inp['mem']['efs'])
    except:
      mem = "500MB"
    mem = std.mem_value(mem)  
    g.memory['efs']['mem'] = mem    
    g.memory['efs']['c'] = int(1 * (mem / 7840))
    g.memory['efs']['g'] = int(12 * (mem / 7840))
    g.memory['efs']['nl'] = int(100 * (mem / 7840))
    
    main.log("EFS " + str(mem))
    main.log("config size       " + str(g.memory['efs']['c']))
    main.log("ghost size        " + str(g.memory['efs']['g']))
    main.log("nl size           " + str(g.memory['efs']['nl']))

    main.log("BP " + str(mem))
    main.log("config size       " + str(g.memory['bp']['c']))
    main.log("ghost size        " + str(g.memory['bp']['g']))
    main.log("nl size           " + str(g.memory['bp']['nl']))

#print(g.memory)

###########################################
#  CLASS label
###########################################
class labels:

  @staticmethod
  def add(label):
    label = label.upper()
    
# Check for mask
    if(label in g.mask.keys()):
      label = g.mask[label].upper()
    
    if(label in g.labels.keys()):
      return label, g.labels[label]
    else: 
      g.labels[label] = len(g.labels) + 1
      return label, g.labels[label]
     
  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'labels.dat', 'w')   
      for k in g.labels.keys():
        fh.write(str(k) + '  ' + str(g.labels[k]) + '\n')
      fh.close()  
        
  @staticmethod
  def get(l_id):     
    for k in g.labels.keys():
      if(l_id == g.labels[k]):
        return k
    return None
        
  @staticmethod
  def add_group(label, group, fn): 
    label = str(label.upper()) 
    group = str(group.upper()) 
    key = label + group
    if(key not in g.grouplabels.keys()):
      g.grouplabels[key] = len(g.grouplabels) + 1
    if(label not in g.groupelement.keys()):
      g.groupelement[label] = [fn]
    else:
      g.groupelement[label].append(fn)
    return key, g.grouplabels[key]
    
  @staticmethod
  def get_group(l_id):     
    for k in g.grouplabels.keys():
      if(l_id == g.grouplabels[k]):
        return k
    return None
        
###########################################
#  CLASS potentia
###########################################
class potential:

  rescale_density_on = False
  rescale_embedding_on = False
  parameters = None

  def run():
  
    efs.init()                           # Initialise (allocate arrays)
    potential.efs_add_potentials()       # Load potentials
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  def load():
    main.log_title("Potential Load")
    
    if(g.inp['potential']['dir'].strip() == ""):
      g.inp['potential']['pot_file'] = g.inp['potential']['index_file']
    else:
      g.inp['potential']['pot_file'] = g.inp['potential']['dir'] + "/" + g.inp['potential']['index_file']
    pot_file = g.inp['potential']['pot_file']
    if(not os.path.isfile(pot_file)):
      main.log("Potential load failed - no pot file")
      return False
    main.log("Potential file: " + str(pot_file))

    if('rescale_density' in g.inp['potential'].keys() and g.inp['potential']['rescale_density'].upper().strip() == "TRUE"):
      potential.rescale_density_on = True

    if('rescale_embedding' in g.inp['potential'].keys() and g.inp['potential']['rescale_embedding'].upper().strip() == "TRUE"):
      potential.rescale_embedding_on = True

# Read potential index
    potential.read_potential(pot_file)
    potential.make_find_fn()
    potential.load_tabulated()
    potential.make_tabulated_points()
    potential.load_analytic()
    potential.make_analytic_points()    
    potential.load_spline()
    potential.make_spline_points()
    potential.load_fit_data()
    potential.rescale_embedding()
    potential.pf_output()
    potential.make_copies()
    potential.parameters = numpy.copy(potential.get_parameters())

    potential.plot_python_potentials(g.dirs['plots'] + "/starting_potential")   
    potential.pf_output_file(g.dirs['input'] + '/start.pot')
#main.end()
       
    return True
  
  @staticmethod
  def make_find_fn():
    for fn in range(len(g.pot_functions['functions'])): 
      t_id = g.pot_functions['functions'][fn]['f_type_id']
      if(g.pot_functions['functions'][fn]['f_type_id'] == 1):
        a_id = g.pot_functions['functions'][fn]['a']
        b_id = g.pot_functions['functions'][fn]['b']
      elif(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        a_id = g.pot_functions['functions'][fn]['a']
        b_id = g.pot_functions['functions'][fn]['f_group']
      elif(g.pot_functions['functions'][fn]['f_type_id'] == 3):
        a_id = g.pot_functions['functions'][fn]['a']
        b_id = g.pot_functions['functions'][fn]['f_group']
      if(t_id not in g.pot_fn.keys()):
        g.pot_fn[t_id] = {}
      if(a_id not in g.pot_fn[t_id].keys()):
        g.pot_fn[t_id][a_id] = {}
      g.pot_fn[t_id][a_id][b_id] = fn
  
  @staticmethod
  def find_fn(t_id, a_id, b_id): 
    if(t_id not in g.pot_fn.keys()):
      return -1
    if(a_id not in g.pot_fn[t_id].keys()):
      return -1
    if(b_id not in g.pot_fn[t_id][a_id].keys()):
      return -1
    return g.pot_fn[t_id][a_id][b_id]

  @staticmethod
  def pot_function():
    return {      
    'f_on': 1,  
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=EMBE, 3=DENS
    'f_group': 1,
    'f_group_label': None,
    'r_cut': 6.5,
    'file': None,
    'function_type': 0,       # 1 tab, 2 analytic
    'f_points': None,         # READ IN TO PYTHON
    'f': None,
    'a_params': None,
    'a_params_fixed': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'a_type': '',
    's_type': '',
    's_nodes': None,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),           # THESE ARE USED BY FORTRAN
    'points_original': numpy.zeros((g.tab_size,g.tab_width,),), 
    'fit_file': None,
    'fit_type': None,         # 1 spline, 2 analytic
    'fit_spline_type': None,  # 1 poly3, 2 poly5
    'fit_parameters': None,
    'fit_parameters_start': None,
    'fit_size': None,
    'fit_mult': None,
    'function_file_name': '',
    }

  @staticmethod
  def rescale_density():
# NOT WORKING YET
    if(potential.rescale_density_on):
      print("Rescale")
      m = {}
      for fn in range(len(g.pot_functions['functions'])): 
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          a_id = g.pot_functions['functions'][fn]['a']
          f_id = g.pot_functions['functions'][fn]['f_group']
          m_rho = rescale_density.estimate_density(fn)
          m[fn] = 1.0 / m_rho
          fn_emb = potential.find_fn(3, a_id, f_id)
          m[fn_emb] = 1.0 / m_rho
          print(m_rho)
# Get parameters
      p = potential.get_parameters()

# Update potential
      a = 0
      for fn in range(len(g.pot_functions['functions'])): 
        multf = 1.0
#if(fn in m.keys()):
#  multf = m[fn]

        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
# Calc b
          g.pot_functions['functions'][fn]['s_nodes'][:,1] = multf * p[a:b]
          potential.make_spline_points_inner(fn)  
        
        elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
# Make Analytic Points
          g.pot_functions['functions'][fn]['a_params'][:] = multf * p[a:b]
          potential.make_analytic_points_inner(fn)

# Update a
        a = b    

# Where the embedding function is analytic, rescale the tab points to range
# from 0.0 up to predicted max density
  @staticmethod
  def rescale_embedding():
    if(potential.rescale_embedding_on):
      m = {}
      for fn in range(len(g.pot_functions['functions'])):         
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):  # Density function
          a_id = g.pot_functions['functions'][fn]['a']
          f_id = g.pot_functions['functions'][fn]['f_group']
          fn_emb = potential.find_fn(3, a_id, f_id)
          m_rho = rescale_density.estimate_density(fn)
          if(fn_emb >= 0):
            m[fn_emb] = m_rho

      for fn in m.keys():
        if(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
          g.pot_functions['functions'][fn]['a_u'] = m[fn_emb]
          potential.make_analytic_points_inner(fn)

  @staticmethod
  def estimate_density(fn):
# Symetric simple cubic used as test case
    sc = {0.86602539999999995: 8, 0.80039053000000004: 24, 0.75: 24, 0.71807032999999998: 24, 0.70710678000000005: 12, 0.72886899000000005: 24, 0.67314560000000001: 48, 0.63737743999999996: 48, 0.625: 24, 0.61237244000000002: 24, 0.57282195999999996: 48, 0.55901699000000005: 24, 0.53033008999999998: 36, 0.51538819999999996: 48, 0.5: 6, 0.64951904999999999: 8, 0.58630196999999995: 24, 0.54486237000000004: 24, 0.46770717000000001: 48, 0.45069390999999998: 24, 0.4145781: 24, 0.39528470999999998: 24, 0.375: 30, 0.43301269999999997: 8, 0.35355339000000002: 12, 0.30618622000000001: 24, 0.27950849999999999: 24, 0.25: 6, 0.21650634999999999: 8, 0.17677670000000001: 12, 0.125: 7}
    a_prim = 0.125
    a = [2.0,3.0,4.0,5.0] 
    m_rho = 0.0
    for a0 in a:
      a0 = (1 / a_prim) * a0
      rho = 0.0    
      for k in sc.keys():
        y = interp.search_x(a0 * k, g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
        rho = rho + sc[k] * y
      m_rho = max(m_rho, rho)
    return m_rho

# LOAD FROM FILE
  @staticmethod
  def read_potential(file_name):
    g.pot_functions['pot_dir'] = ''
    if('/' in file_name):
      lst = file_name.split('/')
      for i in range(len(lst) - 1):
        if(i > 0):
          g.pot_functions['pot_dir'] += '/'
        g.pot_functions['pot_dir'] += lst[i]
    main.log("Loading: " + str(file_name))        

    index = std.config_file_to_list(file_name)  
    pot = potential.pot_function()
    g.pot_functions['pot_file_name'] = file_name
    for row in index:    
      if(len(row) > 1 and row[0].upper() == "POTNAME"):
        g.pot_functions['pot_name'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "ZBLFILE"):
        g.pot_functions['zbl_file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "F_ON"):
        value = row[1].upper()
        if(value[0].upper() == "N"):
          pot['f_on'] = 0
        if(value[0].upper() == "F"):
          pot['f_on'] = 0
      elif(len(row) > 1 and row[0].upper() == "LABEL"):
        label_str, label_id = labels.add(row[1])
        if(label_str not in g.pot_labels):
          g.pot_labels.append(label_str)  # Record potential labels
        pot['a_text'] = label_str
        pot['a'] = label_id
        if(len(row) > 2):
          label_str, label_id = labels.add(row[2])
          if(label_str not in g.pot_labels):
            g.pot_labels.append(label_str)  # Record potential labels
          pot['b_text'] = label_str
          pot['b'] = label_id
      elif(len(row) > 1 and row[0].upper() == "FILE"):
        pot['file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "FIT"):
        pot['fit_file'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "F_TYPE"):
        value = row[1].upper()
        if(value[0] == "P"):
          pot['f_type'] = 'PAIR'
          pot['f_type_id'] = 1
        elif(value[0] == "D"):
          pot['f_type'] = 'DENS'
          pot['f_type_id'] = 2
        elif(value[0] == "E"):
          pot['f_type'] = 'EMBE'
          pot['f_type_id'] = 3
        else:
          pot['f_type'] = 'NONE'
      elif(len(row) > 1 and row[0].upper() == "F_GROUP"):
        pot['f_group_label'] = row[1].strip() 
      elif(len(row) == 1 and row[0].upper() == "F_GROUP"):
        pot['f_group_label'] = "DEFAULT"
      elif(len(row) > 1 and row[0].upper() == "ZOOR"):
        val = row[1].upper() 
        pot['zoor'] = 0  
        if(val[0] == "T" or val[0] == "1" or val[0] == "Y"):
          pot['zoor'] = 1  
      elif(len(row) > 0 and row[0].upper() == "END"):
# END OF POTENTIAL READ
        if(pot['f_type_id'] == 1):
          pot['f_group'] = 0
        else:
          if(pot['a_text'] != None and pot['f_group_label'] != None):
            if(pot['f_group_label'] == ""):
              pot['f_group_label'] = "zzzemptyzzz"
            if(pot['f_group_label'] in g.fgroups.keys()):
              pot['f_group'] = g.fgroups[pot['f_group_label']]
            else:
              id = len(g.fgroups) + 1
              g.fgroups[pot['f_group_label']] = id
              pot['f_group'] = id

        g.pot_functions['functions'].append(pot)
        pot = potential.pot_function()
    
# READ ZBL DATA
    main.log("Load zbl")
    
    read_zbl = False
    for row in index:  
      if(read_zbl == False and row[0].upper() == "ZBLSTART"):  
        read_zbl = True
      elif(read_zbl == True and row[0].upper() == "ZBLEND"):  
        read_zbl = False
      elif(read_zbl):
        label_str, id_1 = labels.add(row[0])
        label_str, id_2 = labels.add(row[1])
        on = True
        if(row[2].upper()[0:1] == "N" or row[2].upper()[0:1] == "F" or row[2].upper()[0:3] == "OFF"):
          on = False
        z1 = float(row[3])
        z2 = float(row[4])
        ra = float(row[5])
        rb = float(row[6])
        
        spline_type = 1
        if(row[7].upper()[0:5] == 'POLY3'):
          spline_type = 1  
        elif(row[7].upper()[0:5] == 'POLY5'):
          spline_type = 2      
        elif(row[7].upper()[0:4] == 'EXP3'):
          spline_type = 3       
        elif(row[7].upper()[0:4] == 'EXP5'):
          spline_type = 4        
        z = {
            'id_1': id_1,
            'id_2': id_2,
            'on': on ,
            'z1': z1 ,
            'z2': z2 ,
            'ra': ra ,
            'rb': rb ,
            'spline_type': spline_type ,
            }
        g.pot_functions['zbl'].append(z)
        main.log("--zbl---")
        main.log(std.dict_to_str(z))

###################
# TABULATED
###################
        
  @staticmethod
  def load_tabulated():      
    for i in range(len(g.pot_functions['functions'])): 
      pf_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][i]['file']      
      if(pf_file is not None and os.path.isfile(pf_file) and potential.pf_file_type(pf_file) == 'T'):
        g.pot_functions['functions'][i]['f_points'] = std.read_csv_array(pf_file, ' ')
        g.pot_functions['functions'][i]['function_type'] = 1
        
# Fill in tabulated points using interp.fill
  @staticmethod
  def make_tabulated_points():  
    for i in range(len(g.pot_functions['functions'])):    
      if(g.pot_functions['functions'][i]['function_type'] == 1):
        g.pot_functions['functions'][i]['points'] = interp.fill(g.pot_functions['functions'][i]['f_points'][:,0], g.pot_functions['functions'][i]['f_points'][:,1], g.tab_size, g.tab_width)   

# Spline and vary accordingly
  @staticmethod
  def vary_tabulated_points(fn, yvar=None):  
    if(type(yvar) != numpy.ndarray):
      yvar = numpy.zeros((10,),)      
    g.pot_functions['functions'][fn]['points'] = spline.vary(g.pot_functions['functions'][fn]['points'][:,0], 
                                                             g.pot_functions['functions'][fn]['points'][:,1], yvar)
                                                             
###################
# SPLINE
###################
  
  @staticmethod
  def load_spline():
    for i in range(len(g.pot_functions['functions'])): 
      pf_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][i]['file']   
      if(pf_file is None):
        print("Error fn " + str(i) + " no file set")
      elif(not os.path.isfile(pf_file)):
        print("Error fn " + str(i) + " file does not exist (" + pf_file + ")")
      elif(pf_file is not None and os.path.isfile(pf_file) and potential.pf_file_type(pf_file) == 'S'):
        g.pot_functions['functions'][i]['function_type'] = 3

        fd = std.config_file_to_list(pf_file)  
#print(fd)
        spline_x = [] 
        spline_y = []
        for l in fd:
          if(l[0].upper() == '#TYPE'):
            g.pot_functions['functions'][i]['s_type'] = l[1].lower()
          elif(l[0].upper().strip() == '#X'):
            for ln in range(1, len(l)):
              p = l[ln].strip()
              if(p != ""):
                try:
                  p = float(p)
                  spline_x.append(p)
                except:
                  pass
          elif(l[0].upper().strip() == '#Y'):
            for ln in range(1, len(l)):
              p = l[ln].strip()
              if(p != ""):
                try:
                  p = float(p)
                  spline_y.append(p)
                except:
                  pass
# Always give a 1 length for param_f even if they aren't used
        g.pot_functions['functions'][i]['s_nodes'] = numpy.zeros((len(spline_x),2,),)

        for p in range(len(spline_x)):
          g.pot_functions['functions'][i]['s_nodes'][p,0] = float(spline_x[p])
          g.pot_functions['functions'][i]['s_nodes'][p,1] = float(spline_y[p])
  
  @staticmethod
  def make_spline_points():  
    for i in range(len(g.pot_functions['functions'])):    
      if(g.pot_functions['functions'][i]['function_type'] == 3):
        potential.make_spline_points_inner(i)
        
  def make_spline_points_inner(fn):  
    st = 1   # poly3
    if(g.pot_functions['functions'][fn]['s_type'] == 'poly5'):
      st = 2
    
    splined = spline.spline_nodes(st, g.pot_functions['functions'][fn]['s_nodes'][:,:], 100)
    g.pot_functions['functions'][fn]['points'] = interp.fill(splined[:,0], splined[:,1], g.tab_size, g.tab_width)
        
# ZBL (if pair)
    potential.make_zbl(fn)

# Now treat as tabulated
    g.pot_functions['functions'][fn]['function_type'] = 1
        
###################
# ANALYTIC
###################
  
  @staticmethod
  def load_analytic():  
    for i in range(len(g.pot_functions['functions'])): 
      pf_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][i]['file']   
      if(pf_file is None):
        print("Error fn " + str(i) + " no file set")
      elif(not os.path.isfile(pf_file)):
        print("Error fn " + str(i) + " file does not exist (" + pf_file + ")")
      elif(pf_file is not None and os.path.isfile(pf_file) and potential.pf_file_type(pf_file) == 'A'):
        g.pot_functions['functions'][i]['function_type'] = 2
        fd = std.config_file_to_list(pf_file)  
#print(fd)
        param = [] 
        param_f = []
        for l in fd:
          if(l[0].upper() == '#TYPE'):
            g.pot_functions['functions'][i]['a_type'] = l[1].lower()
          elif(l[0].upper().strip() == '#P'):
            for ln in range(1, len(l)):
              p = l[ln].strip()
              if(p != ""):
                try:
                  p = float(p)
                  param.append(p)
                except:
                  pass
          elif(l[0].upper().strip() == '#PF'):
            for ln in range(1, len(l)):
              p = l[ln].strip()
              if(p != ""):
                try:
                  p = float(p)
                  param_f.append(p)
                except:
                  pass
          elif(l[0].upper()[0:2] == '#L'):
            g.pot_functions['functions'][i]['a_l'] = float(l[1])
          elif(l[0].upper()[0:2] == '#U'):
            g.pot_functions['functions'][i]['a_u'] = float(l[1])
        g.pot_functions['functions'][i]['a_params'] = numpy.zeros((len(param),),)
        
# Always give a 1 length for param_f even if they aren't used
        if(len(param_f) == 0):
          g.pot_functions['functions'][i]['a_params_fixed'] = numpy.zeros((1,),)
        else:
          g.pot_functions['functions'][i]['a_params_fixed'] = numpy.zeros((len(param_f),),)

        for p in range(len(param)):
          g.pot_functions['functions'][i]['a_params'][p] = float(param[p])
        for p in range(len(param_f)):
          g.pot_functions['functions'][i]['a_params_fixed'][p] = float(param_f[p])

# Save function
        g.pot_functions['functions'][i]['f'] = getattr(potential_functions, g.pot_functions['functions'][i]['a_type'])  
        
  @staticmethod
  def make_analytic_points():  
    for fn in range(len(g.pot_functions['functions'])):  
      if(g.pot_functions['functions'][fn]['function_type'] == 2):  
        potential.make_analytic_points_inner(fn)
        
  @staticmethod
  def make_analytic_points_inner(fn):  
# Temp x,y array
    temp = numpy.zeros((g.tab_size,2,),)
    temp[:,0] = numpy.linspace(g.pot_functions['functions'][fn]['a_l'], g.pot_functions['functions'][fn]['a_u'], g.tab_size)
    temp[:,1] = g.pot_functions['functions'][fn]['f'](temp[:,0],  
                                                      g.pot_functions['functions'][fn]['a_params'],  
                                                      g.pot_functions['functions'][fn]['a_params_fixed'])
       
# Interpfill
    g.pot_functions['functions'][fn]['points'] = interp.fill(temp[:,0],temp[:,1], g.tab_size, g.tab_width)
    
# ZBL (if pair)
    potential.make_zbl(fn)
    
  @staticmethod
  def load_fit_data(): 
# READ FIT DATA
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_file'] != None):
        params = None
        fit_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][fn]['fit_file']
        if(os.path.isfile(fit_file)):
# Read File
          params_start = None
          fh = open(fit_file, 'r')
          for line in fh:
            line = std.one_space(line.strip().upper())
            f = line.split(" ")
            if(line[0:4] == "#FIT"):
              if(f[-1] == 'S'):
                g.pot_functions['functions'][fn]['fit_type'] = 1
              elif(f[-1] == 'A'):
                g.pot_functions['functions'][fn]['fit_type'] = 2                
            if(line[0:3] == "#PL"):
              params_lower = f[1:]                     
            if(line[0:3] == "#PU"):
              params_upper = f[1:]                     
            if(line[0:3] == "#PS"):
              params_start = f[1:]                      
            if(line[0:2] == "#N"):
              spline_nodes = f[1:]                      
            if(line[0:3] == "#ST"):
              spline_type = 1
              if(f[1].upper().strip() == "POLY5"):
                spline_type = 2
              
          fh.close()
      
          if(g.pot_functions['functions'][fn]['fit_type'] != None and type(params_lower) == list):
          
# g.pot_functions['functions'][fn]['fit_type']       1 = Spline  2 = Analytic
# g.pot_functions['functions'][fn]['function_type']  1 = Tab     2 = Analytic
        
# TAB DATA BUT WRONG FITTING
            if(g.pot_functions['functions'][fn]['fit_type'] == 2 and
               g.pot_functions['functions'][fn]['function_type'] == 1):
              g.pot_functions['functions'][fn]['fit_type'] = 1
              g.pot_functions['functions'][fn]['fit_size'] = 15
              g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((3,15,),)  
              for i in range(g.pot_functions['functions'][fn]['fit_size']):
                g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(-1.0)
                g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(1.0)
                
# Analytic fitting
            elif(g.pot_functions['functions'][fn]['fit_type'] == 2 and
               g.pot_functions['functions'][fn]['function_type'] == 2):
              g.pot_functions['functions'][fn]['fit_size'] = len(g.pot_functions['functions'][fn]['a_params'])
              g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((3,len(g.pot_functions['functions'][fn]['a_params']),),) 
              g.pot_functions['functions'][fn]['fit_parameters_start'] = numpy.zeros((len(g.pot_functions['functions'][fn]['a_params']),),)    
            
              for i in range(g.pot_functions['functions'][fn]['fit_size']):
                g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
                g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])
#g.pot_functions['functions'][fn]['fit_parameters_start'][i] = float(params_start[i])
                
                if(params_start is None):
                  g.pot_functions['functions'][fn]['fit_parameters_start'][i] = 0.5 * (float(params_lower[i]) + float(params_upper[i]))
                else:
                  g.pot_functions['functions'][fn]['fit_parameters_start'][i] = float(params_start[i])
              
# Spline fitting
            elif(g.pot_functions['functions'][fn]['fit_type'] == 1):
              g.pot_functions['functions'][fn]['fit_size'] = len(params_lower)
              g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((3,len(params_lower),),)
              g.pot_functions['functions'][fn]['fit_spline_type'] = spline_type  
              for i in range(g.pot_functions['functions'][fn]['fit_size']):
                g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
                g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])
                g.pot_functions['functions'][fn]['fit_parameters'][2,i] = float(spline_nodes[i])

  @staticmethod
  def make_copies(): 
    for fn in range(len(g.pot_functions['functions'])): 
      g.pot_functions['functions'][fn]['points_original'] = numpy.copy(g.pot_functions['functions'][fn]['points'])
  
    g.pot_functions['functions_original'] = copy.deepcopy(g.pot_functions['functions'])

  @staticmethod
  def spline_prep(): 
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 2):   # If spline fit
        yvar = numpy.zeros((len(g.pot_functions['functions'][fn]['fit_parameters']),),)
        potential.vary_tabulated_points(fn, yvar)
    
  @staticmethod
  def pf_file_type(file): 
    fh = open(file, 'r')
    for row in fh:
      fh.close()
      if(row.strip().upper()[0:2] == '#A'):
        return 'A'
      elif(row.strip().upper()[0:2] == '#S'):
        return 'S'
      return 'T'
      
  @staticmethod
  def pf_output(): 
    if(g.outputs):     
      potential.pf_output_file(g.dirs['output'] + '/' + 'pot.dat')
    
  @staticmethod
  def pf_output_file(file_path): 

      fh = open(file_path, 'w')
      
      n = 0
      for pf in g.pot_functions['functions']:
        n = n + 1
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('POTENTIAL ' + str(n) + '\n')
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('f_on: ' + str(pf['f_on']) + '\n')
        fh.write('a_text: ' + str(pf['a_text']) + '\n')
        fh.write('b_text: ' + str(pf['b_text']) + '\n')
        fh.write('a: ' + str(pf['a']) + '\n')
        fh.write('b: ' + str(pf['b']) + '\n')
        fh.write('f_type: ' + str(pf['f_type']) + '\n')
        fh.write('f_group: ' + str(pf['f_group']) + '\n')
        fh.write('f_group_label: ' + str(pf['f_group_label']) + '\n')
        fh.write('r_cut: ' + str(pf['r_cut']) + '\n')
        fh.write('file: ' + str(pf['file']) + '\n')
        fh.write('function_type: ' + str(pf['function_type']) + '\n')
        fh.write('a_type: ' + str(pf['a_type']) + '\n')
        fh.write('f: ' + str(pf['f']) + '\n')
        fh.write('a_l: ' + str(pf['a_l']) + '\n')
        fh.write('a_u: ' + str(pf['a_u']) + '\n')
        fh.write('zoor: ' + str(pf['zoor']) + '\n')
        fh.write('a_params: \n')
        try:
          for i in range(len(pf['a_params'])):
            fh.write(str(pf['a_params'][i]) + '\n')    
        except:
          fh.write(str('None\n')) 
          
        fh.write('f_points: \n')
        try:
          for i in range(len(pf['f_points'])):
            for j in range(len(pf['f_points'][i])):
              fh.write(str(pf['f_points'][i, j]) + ' ')  
            fh.write('\n')  
        except:
          fh.write(str('None\n')) 
          
        fh.write('points: \n')
        try:
          for i in range(len(pf['points'])):
            for j in range(len(pf['points'][i])):
              fh.write(str(pf['points'][i, j]) + ' ')  
            fh.write('\n') 
        except:
          fh.write(str('None\n'))         
      fh.close()
      
  @staticmethod
  def get_parameters():  
    w = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_size'] != None):
        w = w + g.pot_functions['functions'][fn]['fit_size']
    p = numpy.zeros((w,),)  
      
# Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_size'] != None):
        b = a + g.pot_functions['functions'][fn]['fit_size']  
        if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
          p[a:b] = g.pot_functions['functions'][fn]['s_nodes'][:,1]
        elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
          p[a:b] = g.pot_functions['functions'][fn]['a_params'][:]
# Update a
        a = b    
    
    return p
    
  @staticmethod
  def get_start_parameters():   
    w = 0
    for fn in range(len(g.pot_functions['functions'])): 
      w = w + g.pot_functions['functions'][fn]['fit_size']
    p = numpy.zeros((w,),)  
      
# Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      b = a + g.pot_functions['functions'][fn]['fit_size']  
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
        p[a:b] = g.pot_functions['functions'][fn]['fit_parameters_start'][:,1]
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        p[a:b] = g.pot_functions['functions'][fn]['fit_parameters_start'][:]
# Update a
      a = b    
      
    return p
    
  @staticmethod
  def save_potential(dir_out=None, name_out=None):
    potential_output.eampa(dir_out, name_out)
      
  """
  def pot_function():
    return {      
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
    }
  """    
      
  def plot_python_potentials(dir = None):  
  
    try:
      if(dir == None):
        dir = g.dirs['pots']
      std.make_dir(dir)
      print(dir)
      
      for i in range(len(g.pot_functions['functions'])): 
        pot_name = 'py_pot'
        pot_count = i + 1
      
        pot_type = g.pot_functions['functions'][i]['f_type_id']
        label_a = g.pot_functions['functions'][i]['a']
      
        if(pot_type == 1):
          label_b = g.pot_functions['functions'][i]['b']
        else:
          label_b = g.pot_functions['functions'][i]['f_group']
        
        potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                                 g.pot_functions['functions'][i]['points'][:,0], 
                                 g.pot_functions['functions'][i]['points'][:,1], 
                                 g.pot_functions['functions'][i]['points'][:,2], 
                                 g.pot_functions['functions'][i]['points'][:,3])
      return True
    except:
      return False  
      
  def plot_fortran_potentials(dir = None):  
    
    try:
    
      if(dir == None):
        dir = g.dirs['pots']
      std.make_dir(dir)
      print(dir)
    
      if(efs.pc > 0):
        for i in range(efs.pc):
    
          a = efs.pkey[i,0] - 1
          b = efs.pkey[i,1]
      
          pot_name = 'fort_pot'
          pot_count = i + 1
       
          pot_type = efs.pkey[i,2]
          label_a = efs.pkey[i,3]
          label_b = efs.pkey[i,4]
      
          potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                                   efs.pot[a:b,0], efs.pot[a:b,1], 
                                   efs.pot[a:b,2], efs.pot[a:b,3])
    
      elif(bp.pc > 0):
        for i in range(bp.pc):
    
          a = bp.pkey[i,0] - 1
          b = bp.pkey[i,1]
      
          pot_name = 'fort_pot'
          pot_count = i + 1
       
          pot_type = bp.pkey[i,2]
          label_a = bp.pkey[i,3]
          label_b = bp.pkey[i,4]
      
          potential.plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, 
                                   bp.pot[a:b,0], bp.pot[a:b,1], 
                                   bp.pot[a:b,2], bp.pot[a:b,3])
      return True
    except:
      return False                                  
                                 
  def plot_potential(dir, pot_name, pot_count, pot_type, label_a, label_b, x, y, yp, ypp):  
      
      if(pot_type == 1):
        plot_title = "Pair - " + labels.get(label_a) + ' ' + labels.get(label_b)    
        x_axis = "Seperation (Ang)"
        y_axis = "Potential Energy (eV)"
      elif(pot_type == 2):   
        plot_title = "Density - " + labels.get(label_a) + ' Group ' + str(label_b)    
        x_axis = "Seperation (Ang)"
        y_axis = "Electron Density"   
      elif(pot_type == 3):
        plot_title = "Embedding - " + labels.get(label_a) + ' Group ' + str(label_b)     
        x_axis = "Electron Density"
        y_axis = "Potential Energy (eV)"
      
      file_name = str(pot_count)
      while(len(file_name)<3):
        file_name = '0' + file_name
      file_name = pot_name + '_' + file_name + '.eps'
      
      plt.clf()    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')
      fig, axs = plt.subplots(1, 1, figsize=(12,9))
      fig.tight_layout(pad=5.0)   
      
      fig.suptitle(plot_title)
      axs.plot(x, y, color='k', ls='solid')
      axs.plot(x, yp, color='k', ls='dashed')
      axs.plot(x, ypp, color='k', ls='dotted')
      axs.set_title('')
      axs.set_xlabel(x_axis)
      axs.set_ylabel(y_axis)
      
      min_y = min(y)
      max_y = max(y)
      r = max_y - min_y
      
      min_y = min_y - 0.05 * r
      max_y = max_y + 0.05 * r
      
      if(min_y < -1.0E03):
        min_y = -1.0E03
      if(max_y > 1.0E04):
        max_y = 1.0E04
         
      axs.set_ylim(min_y, max_y)
      axs.set_yscale('symlog', linthreshy=10)
      plt.savefig(dir + '/' + file_name, format='eps')
      plt.close('all') 
      
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
      
  def parameter_count():
# Find width
    p_count = 0
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        p_count = p_count + g.pot_functions['functions'][fn]['fit_size'] 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        p_count = p_count + g.pot_functions['functions'][fn]['fit_size']
    return p_count  
      
###########################################################
# ZBL
###########################################################
      
  def make_zbl(fn):
    
    if(g.pot_functions['functions'][fn]['f_type_id'] == 1):
      zbl_n = -1
      if(len(g.pot_functions['zbl'])>0):
        for i in range(len(g.pot_functions['zbl'])): 
          if((g.pot_functions['functions'][fn]['a'] == g.pot_functions['zbl'][i]['id_1'] 
             and g.pot_functions['functions'][fn]['b'] == g.pot_functions['zbl'][i]['id_2'])
             or (g.pot_functions['functions'][fn]['a'] == g.pot_functions['zbl'][i]['id_2'] 
             and g.pot_functions['functions'][fn]['b'] == g.pot_functions['zbl'][i]['id_1'])):
            zbl_n = i
            break
      
      if(zbl_n>=0):
        ra = g.pot_functions['zbl'][zbl_n]['ra']
        rb = g.pot_functions['zbl'][zbl_n]['rb']
        qa = g.pot_functions['zbl'][zbl_n]['z1']
        qb = g.pot_functions['zbl'][zbl_n]['z2']
        spline_type = g.pot_functions['zbl'][zbl_n]['spline_type']
        
# Find the point that need replacing
        ra_n = -1
        for n in range(len(g.pot_functions['functions'][fn]['points'])):
          if(g.pot_functions['functions'][fn]['points'][n,0] <= ra):
            ra_n = n
          elif(g.pot_functions['functions'][fn]['points'][n,0] <= rb):
            rb_n = n
            
# Get node values
        ya = fnc.zbl(ra, [qa, qb], [0.0])  
        yap = fnc.zbl_dydr(ra, [qa, qb], [0.0]) 
        yb = interp.interpolate(rb, g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1], 4)
        ybp = interp.interpolate_dydxn(rb, g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1], 1, 4)

# Get spline coeffs
        coeffs = spline.spline_ab(spline_type, [ra, ya, yap], [rb, yb, ybp])

# Replace with ZBL and SPLINE
        g.pot_functions['functions'][fn]['points'][:ra_n,1] = fnc.zbl_v(g.pot_functions['functions'][fn]['points'][:ra_n,0], [qa, qb], [0.0])        
        g.pot_functions['functions'][fn]['points'][ra_n:rb_n,1] = potential.poly(g.pot_functions['functions'][fn]['points'][ra_n:rb_n,0], coeffs)
        
# Recalculate derivatives etc
        g.pot_functions['functions'][fn]['points'] = interp.fill(g.pot_functions['functions'][fn]['points'][:,0],g.pot_functions['functions'][fn]['points'][:,1], g.tab_size, g.tab_width)
    
  def poly(x, c):
    y = 0.0
    for i in range(len(c)):
      y = y + c[i] * x**i
    return y
      
###########################################################
# F2PY functions
###########################################################

  def efs_add_potentials():
  
# Clear
    efs.clear_potentials()
    
# ADD POTENTIALS
    for pn in range(len(g.pot_functions['functions'])):
      pf = g.pot_functions['functions'][pn]
#print(pf['f_type_id'], pf['a'], pf['f_group'])
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          fortran_id = efs.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] ,
                            pn
                           )
          g.pot_functions['functions'][pn]['fortran_id'] = fortran_id
        elif(pf['f_type_id'] > 1):
          fortran_id = efs.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'],
                            pn
                           )   
          g.pot_functions['functions'][pn]['fortran_id'] = fortran_id
    
# SET POTENTIALS
    efs.set_potentials() 
#exit()
    
  def es_add_potentials():
  
# Clear
    es.clear_potentials()
    
# ADD POTENTIALS
    for pn in range(len(g.pot_functions['functions'])):
      pf = g.pot_functions['functions'][pn]
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          es.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          es.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )    
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      es.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    

# SET POTENTIALS
    es.set_potentials()     

  def bp_add_potentials():
  
# Clear
    bp.clear_potentials()
    
# ADD POTENTIALS
    for pn in range(len(g.pot_functions['functions'])):
      pf = g.pot_functions['functions'][pn]
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          fortran_id = bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'],
                            pn
                           )
        elif(pf['f_type_id'] > 1):
          fortran_id = bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'],
                            pn 
                           )    
      
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      bp.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    
                                               
# SET POTENTIALS
    bp.set_potentials()
      
  def relax_add_potentials():
  
# Clear
    relax.clear_potentials()
    
# ADD POTENTIALS
    for pf in g.pot_functions['functions']:
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          relax.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          relax.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )    
      
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      relax.add_zbl(
                  zbl['id_1'],
                  zbl['id_2'], 
                  zbl['on'],  
                  zbl['z1'], 
                  zbl['z2'] , 
                  zbl['ra'] , 
                  zbl['rb'] , 
                  zbl['spline_type'] 
                 )    
                                               
# SET POTENTIALS
    relax.set_potentials()
      
###########################################
#  CLASS potential_function
###########################################
class potential_functions:

##############################
#  All the functions have the same format
#  r which will be a vector/1D array of values
#  p which will be a vector of parameters
#  pf which will be a vector of "fixed" parameters - such as r_cut which will be set and fixed, and not varied during the fitting process

##############################################
# PAIR FUNCTIONS
##############################################

# Lennard Jones Potential
# p[0] = e
# p[1] = rm
# f(x) = A * ((B / r)**12 - 2 * (B/r)**6)
  @staticmethod
  def lennard_jones(r, p, pf):
    return fnc.lennard_jones_v(r, p, pf)
    
# Morse Potential
# p[0] = d
# p[1] = a
# p[2] = re
# f(x) = A * (exp(-2.0D0 * B * (r - C)) - 2.0D0 * exp(-B*(r - C)))
  @staticmethod
  def morse(r, p, pf):
    return fnc.morse_v(r, p, pf)
#return

# Buckingham Potential
# p[0] = A
# p[1] = B
# p[2] = C
# f(x) = A * exp(-1 * B * r) - C / r**6
  @staticmethod
  def buckingham(r, p, pf):
    return fnc.buckingham_v(r, p, pf)

  @staticmethod
  def ackland_mendelev_pair(r, p, pf):
    return fnc.ackland_mendelev_pair_v(r, p, pf)

  @staticmethod
  def pair_spline(r, p, pf):
    return fnc.pair_spline_v(r, p, pf)

##############################################
# DENSITY FUNCTIONS
##############################################

# Embedding Finnis-Sinclair
  @staticmethod
  def quadratic_density(r, p, pf):
    return fnc.quadratic_density_v(r, p, pf)
    
##############################################
# EMBEDDING FUNCTIONS
##############################################

# Embedding Finnis-Sinclair
# f(x) = -A sqrt(rho)
  @staticmethod
  def fs_embedding(r, p, pf):
    return fnc.fs_embedding(r, p, pf)

# Embedding Mendelev
# f(x) = -sqrt(rho) + A*rho**2
  @staticmethod
  def mendelev_embedding(r, p, pf):
    return fnc.mendelev_embedding_v(r, p, pf)

# Triple Embedding
# f(x) = A * sqrt(r) + B * r + C * r**2
  @staticmethod
  def triple_embedding(r, p, pf):
    return fnc.triple_embedding_v(r, p, pf)

# Triple Embedding
# f(x) = A + B * sqrt(r) + C * r**2 + D * r**4
  @staticmethod
  def quad_embedding(r, p, pf):
    return fnc.quad_embedding_v(r, p, pf)
    
# Embedding Ackland (Olsson/Walenius)
# f(x) = A sqrt(rho) + B rho**2 + C rho**4
  @staticmethod
  def ackland_embedding(r, p, pf):
    return fnc.ackland_embedding_v(r, p, pf)

##############################################
# SIMPLE SPLINES
##############################################

  @staticmethod
  def cubic_spline(r, p, pf):
    return fnc.cubic_spline_v(r, p, pf)
    
  @staticmethod
  def quintic_spline(r, p, pf):
    return fnc.quintic_spline_v(r, p, pf)
    
##############################################
# NODE SPLINES
##############################################
   
  @staticmethod
  def spline_n_node(r, p, pf):
    return fnc.spline_n_node_v(r, p)
    
  @staticmethod
  def cubic_knot_spline(r, p, pf):
    return fnc.cubic_knot_spline_v(r, p, pf)
    
  @staticmethod
  def cubic_knot_spline_fixed_end(r, p, pf):
    return fnc.cubic_knot_spline_fixed_end_v(r, p, pf)

###########################################
#  CLASS potential_outpu
###########################################
class potential_output:

  @staticmethod
  def full(dir_out = None):
  
    if(dir_out==None):  
      dir_out = g.dirs['wd'] + '/' + g.fit_results['results_dir'] 
  
# Make output directory
    std.make_dir(dir_out)
    
    potential_output.best_parameters(dir_out) 
    potential_output.eampa(dir_out)
    potential_output.data_file(dir_out)
    potential_output.dl_poly(dir_out)
    potential_output.plots(dir_out)
    
  def best_parameters(dir_out = None):
    if(dir_out==None):  
      dir_out = g.dirs['results'] + '/best_parameters'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/best_parameters.txt', 'w')    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        fh.write("Fn: " + str(fn) + "[S]  ")
        fh.write("\n")
        for i in range(len(g.pot_functions['functions'][fn]['s_nodes'][:,1])):
          fh.write("  [" + str(i) + "] " + display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,0],8) + " ")
        for i in range(len(g.pot_functions['functions'][fn]['s_nodes'][:,1])):
          fh.write("  [" + str(i) + "] " + display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,1],8) + " ")
          
        fh.write("\n")
#for i in range(g.pot_functions['functions'][fn]['fit_size']):
#  fh.write("[" + str(i) + "]" + display.pad_r(g.pot_functions['functions'][fn]['fit_parameters'][0,i],8))
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        fh.write("Fn:" + str(fn) + "[A]  ")
        fh.write("\n")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          fh.write("  [" + str(i) + "] " + display.pad_r(g.pot_functions['functions'][fn]['a_params'][i],8) + " ")
        fh.write("\n")
    fh.close()
       
  @staticmethod
  def eampa(dir_out = None):
    if(dir_out==None):   
      dir_out = g.dirs['results'] + '/pot_save'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/out.pot', 'w')
    fh.write('POTNAME ' + '\n')  
    fh.write('\n')  
    fh.write('\n')  
    
    for fn in range(len(g.pot_functions['functions'])): 
      f = g.pot_functions['functions'][fn]

# File Name
#pf_name = f['f_type'] + '_' + f['a_text'].strip()
#if(f['b_text'].strip() != ''):
#  pf_name += '_' + f['b_text'].strip()
#if(str(f['f_group']).strip() != ''):
#  pf_name += '_' + str(f['f_group'])
#pf_name = pf_name.lower() + '.pot'
      pf_name = f['file']
      
      a = labels.get(f['a'])
      b = labels.get(f['b'])
      fh.write('! ' + f['f_type'] + '  ' + f['a_text'] + '  ' + f['b_text'] + '\n')
      fh.write('START\n')
      if(f['f_on'] == 1):
        fh.write('F_ON true\n')
      else:
        fh.write('F_ON false\n')
      fh.write('F_TYPE ' + f['f_type'] + '\n')
      fh.write('FILE ' + pf_name + '\n')
      fh.write('LABEL ' + f['a_text'] + '  ' + f['b_text'] + '\n')
      fh.write('F_GROUP ' + str(f['f_group']) + '\n')
      fh.write('R_CUT ' + str(f['r_cut']) + '\n')
      fh.write('END\n')  
      fh.write('\n')  
      fh.write('\n')  
      
# Write Function File
      fh_pf = open(dir_out + '/' + pf_name, 'w')
      if(f['function_type'] == 1):
        pass
      elif(f['function_type'] == 2):
        fh_pf.write('#A\n')
        fh_pf.write('#TYPE ' + str(f['a_type']) + '\n')
        fh_pf.write('#P')
        if(f['a_params'] is not None):
          for p in f['a_params']:        
            fh_pf.write(' ' + str(p))
          fh_pf.write('\n')
        if(f['a_params_fixed'] is not None):
          fh_pf.write('#PF')
          for p in f['a_params_fixed']:        
            fh_pf.write(' ' + str(p))
          fh_pf.write('\n')
        fh_pf.write('#L ' + str(f['a_l']) + '\n')
        fh_pf.write('#U ' + str(f['a_u']) + '\n')
      fh_pf.close()

# Write Function File
      if(f['function_type'] == 2):
        fh_pfit = open(dir_out + '/' + f['fit_file'], 'w')
        fh_pfit.write('#FIT A' + '\n')
        fh_pfit.write('#PS ')
        for p in f['a_params']:        
          fh_pfit.write(' ' + str(p))
        fh_pfit.write('\n')
        fh_pfit.write('#PL ')
        for i in range(len(f['a_params'])):    
          fh_pfit.write(' ' + str(str(f['a_params'][i] - 0.5 * (f['fit_parameters'][1,i] - f['fit_parameters'][0,i]))))
        fh_pfit.write('\n')
        fh_pfit.write('#PU ')
        for i in range(len(f['a_params'])):    
          fh_pfit.write(' ' + str(str(f['a_params'][i] + 0.5 * (f['fit_parameters'][1,i] - f['fit_parameters'][0,i]))))
        fh_pfit.write('\n')
        fh_pfit.close()
        
    fh.close()
  
  @staticmethod
  def data_file(dir_out = None):
    if(dir_out==None):  
      dir_out = g.dirs['results'] + '/pot_fortran'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/data.pot', 'w')
    if(efs.pc > 0):
      for i in range(efs.pc):    
        a = efs.pkey[i,0] - 1
        b = efs.pkey[i,1]       
        pot_type = efs.pkey[i,2]
        label_a = efs.pkey[i,3]
        label_b = efs.pkey[i,4]        
        fh.write(str(pot_type) + " " + str(label_a) + " " + str(label_b) + "\n")
        for n in range(a,b,1):
          fh.write(potential_output.pad_r(efs.pot[n,0], 20) + " ")
          fh.write(potential_output.pad_r(efs.pot[n,1], 20) + " ")
          fh.write(potential_output.pad_r(efs.pot[n,2], 20) + " ")
          fh.write(potential_output.pad_r(efs.pot[n,3], 20) + "\n")  
    fh.close()
    
  @staticmethod
  def dl_poly(dir_out = None):
    if(dir_out==None):  
      dir_out = g.dirs['results'] + '/dl_poly'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/pot.eam', 'w')
    fh.write('# EAMPA Output DL_POLY\n')
    if(efs.pc > 0):
      fh.write(str(efs.pc) + '\n')
      for i in range(efs.pc):    
        a = efs.pkey[i,0] - 1
        b = efs.pkey[i,1]       
        pot_type = efs.pkey[i,2]
        label_a = labels.get(efs.pkey[i,3])
        label_b = labels.get(efs.pkey[i,4]) 
        pykey = efs.pkey_python[i]
        pf = g.pot_functions['functions'][pykey]
        if(pot_type == 1):
          fh.write('pair ')   
          fh.write(label_a + ' ')   
          fh.write(label_b + ' ')  
        elif(pot_type == 2):
          if(pf['f_group_label'] == "1"):
            fh.write('dden ')  
          elif(pf['f_group_label'] == "2"):
            fh.write('sden ')  
          else:
            fh.write('dens ')  
          fh.write(label_a + ' ')  
        elif(pot_type == 3):
          if(pf['f_group_label'] == "1"):
            fh.write('demb ')  
          elif(pf['f_group_label'] == "2"):
            fh.write('semb ')  
          else:
            fh.write('embe ')  
          fh.write(label_a + ' ')   
        fh.write(str(b - a) + ' ') 
        fh.write(str(efs.pot[a,0]) + ' ') 
        fh.write(str(efs.pot[b-1,0]) + ' ') 
        fh.write('\n')
        
        ended = False
        for n in range(a,b,1):
          ended = False
          fh.write(potential_output.pad_r("{:.14e}".format(efs.pot[n,1]), 20) + " ")
          if(n % 4 == 3):
            fh.write('\n')
            ended = True
            
        if(not ended):
          fh.write('\n')
          
    fh.close()
    
  @staticmethod
  def plots(dir_out = None): 
    if(dir_out==None):  
      dir_out = g.dirs['results'] + '/plots'
    std.make_dir(dir_out)    
    potential.plot_python_potentials(dir_out)
    potential.plot_fortran_potentials(dir_out)
  
  @staticmethod
  def pad_r(inp, p=7):
    if(inp == None):
      return ""    
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]

###########################################
#  CLASS rescale_densit
###########################################
class rescale_density:
  
  def run():
    m = {}
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        a_id = g.pot_functions['functions'][fn]['a']
        f_id = g.pot_functions['functions'][fn]['f_group']
        m_rho = rescale_density.estimate_density(fn)
        m[fn] = 1.0 / m_rho
        fn_emb = potential.find_fn(3, a_id, f_id)
        m[fn_emb] = 1.0 / m_rho

    """
# Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
# Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size']  
        
        g.pot_functions['functions'][fn]['s_nodes'][:,1] = p[a:b]
        potential.make_spline_points_inner(fn)
        
# Make Analytic Points
#g.pot_functions['functions'][fn]['s_params'][:] = p[a:b]
# LOAD ORIGINAL
#g.pot_functions['functions'][fn]['points'] = numpy.copy(g.pot_functions['functions'][fn]['points_original'])
# VARY SPLINE
#potential.vary_tabulated_points(fn, p[a:b])
        
# Update a
        a = b   
        
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
# Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
# Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)
        a = b    
    """

  def estimate_densities():
    out = []
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        out.append([fn, rescale_density.estimate_density(fn)])
    return out

  def estimate_density(fn):
# Symetric simple cubic used as test case
    sc = {0.86602539999999995: 8, 0.80039053000000004: 24, 0.75: 24, 0.71807032999999998: 24, 0.70710678000000005: 12, 0.72886899000000005: 24, 0.67314560000000001: 48, 0.63737743999999996: 48, 0.625: 24, 0.61237244000000002: 24, 0.57282195999999996: 48, 0.55901699000000005: 24, 0.53033008999999998: 36, 0.51538819999999996: 48, 0.5: 6, 0.64951904999999999: 8, 0.58630196999999995: 24, 0.54486237000000004: 24, 0.46770717000000001: 48, 0.45069390999999998: 24, 0.4145781: 24, 0.39528470999999998: 24, 0.375: 30, 0.43301269999999997: 8, 0.35355339000000002: 12, 0.30618622000000001: 24, 0.27950849999999999: 24, 0.25: 6, 0.21650634999999999: 8, 0.17677670000000001: 12, 0.125: 7}
    a_prim = 0.125
    a = [2.0,3.0,4.0,5.0] 
    m_rho = 0.0
    for a0 in a:
      a0 = (1 / a_prim) * a0
      rho = 0.0    
      for k in sc.keys():
        y = interp.search_x(a0 * k, g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
        rho = rho + sc[k] * y
      m_rho = max(m_rho, rho)
    return m_rho
    
  def max_densities():    
    density_list = None
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        rho = rescale_density.estimate_density(fn) 
        if(density_list is None):
          density_list = str(rho)
        else:
          density_list = density_list + "/" + str(rho)
    return density_list
        
###########################################
#  CLASS potential_var
###########################################
class potential_vary:
 
  def vary_all():
    for pn in range(len(g.pot_functions['functions'])):
      if(g.pot_functions['functions'][pn]['function_type'] == 1):
        potential_vary.vary_tabulated(pn)
      elif(g.pot_functions['functions'][pn]['function_type'] == 2):
        potential_vary.vary_analytic(pn)
  
  def vary_tabulated(pn):
    pass
  
  def vary_analytic(pn):
    print(g.pot_functions['functions'][pn])
    
###########################################
#  CLASS config
###########################################
class configs:

  @staticmethod
  def load():
    try:
      configs_dir = g.inp['configs']['dir']
    except:  
      return False
    if(not os.path.isdir(configs_dir)):
      return False      
    
# Make file list
    g.configs['config_files'] = []
    g.configs['config_files'] = configs.load_files(configs_dir, g.configs['config_files'])
    
# Read files in
    configs.read()

    return True

  @staticmethod
  def load_files(path_in, files):
    for path in os.listdir(path_in):
      path_new = path_in + "/" + path
      if(os.path.isdir(path_new)):
        files = configs.load_files(path_new, files)
      else:
        type = configs.file_type(path_new)
        files.append([path_new,type])
    return files
    
  @staticmethod
  def read():
# Loop through files and read in
    for file in g.configs['config_files']: 
      main.log(file[1] + "  " + file[0])
      if(file[1] == 'std'):        
        configs.std(file[0])
      elif(file[1] == 'qe'):
        configs.qe(file[0])
      elif(file[1] == 'raw'):
        configs.raw(file[0])
        
  @staticmethod
  def std(file_path):   
  
# Read content from file
    content = std.file_to_list(file_path)
    content = std.prep_data(content)
    content = std.remove_comments(content)
  
# Split up into individual configs
    config_list = []    
    temp = []
    last_line = None
    for line in content:
      if(last_line != None and line != '' and last_line != '' and last_line[0] != "#" and line[0] == "#"):
        config_list.append(temp)
        temp = []
      if(line != ''):  
        temp.append(line)      
        last_line = line
    if(len(temp)>0):
      config_list.append(temp)
      
# Read each config
    for i in range(len(config_list)):
      configs.add_config(config_list[i], file_path, i)

  @staticmethod
  def make_config():
    return {
    'f_id': None, 
    'file_type': '',  
    'file_path': '',  
    'file_part': 0,  
    'alat': 0.0,  
    'uv_prim': numpy.zeros((3,3,),),
    'uv': numpy.zeros((3,3,),),
    'stress': numpy.zeros((3,3,),),
    'c': numpy.zeros((3,),dtype=numpy.int32,),
    'rcut': 0.0,
    'rverlet': 0.0,
    'mtemp': {},
    'm': {},
    'coord_count_prim': 0,
    'coords_label_prim': None,
    'coords_label_id_prim': None,
    'coords_prim': None,
    'forces_prim': None,
    'coord_count': 0,
    'coords_label': None,
    'coords_label_id': None,
    'coords': None,
    'forces': None,
    'energy_per_atom': None,
    'energy': None,
    'e': 0,
    'f': 0,
    's': 0,
    'n_atoms_prim': 0,
    'n_atoms': 0,
    'l_units': 'ANG',
    'e_units': 'EV',
    'f_units': 'EV/ANG',
    's_units': 'EV/ANG3',
    }    
  
  @staticmethod
  def add_config(content, file_path, i, config_type='STANDARD'):
    fd = configs.make_config()
    fd['file_type'] = config_type
    fd['file_path'] = file_path
    fd['file_part'] = i

    fd['c'][0] = 1
    fd['c'][1] = 1
    fd['c'][2] = 1
    
# FIRST READ
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      epa_set = False

      if(len(f) > 1 and f[0].upper() == "#ALAT"):
        fd['alat'] = float(f[1])
      if(len(f) > 3 and f[0].upper() == "#X"):
        fd['uv_prim'][0,0] = float(f[1])
        fd['uv_prim'][0,1] = float(f[2])
        fd['uv_prim'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Y"):
        fd['uv_prim'][1,0] = float(f[1])
        fd['uv_prim'][1,1] = float(f[2])
        fd['uv_prim'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#Z"):
        fd['uv_prim'][2,0] = float(f[1])
        fd['uv_prim'][2,1] = float(f[2])
        fd['uv_prim'][2,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SX"):
        fd['s'] = 1
        fd['stress'][0,0] = float(f[1])
        fd['stress'][0,1] = float(f[2])
        fd['stress'][0,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SY"):
        fd['s'] = 1
        fd['stress'][1,0] = float(f[1])
        fd['stress'][1,1] = float(f[2])
        fd['stress'][1,2] = float(f[3])
      if(len(f) > 3 and f[0].upper() == "#SZ"):
        fd['s'] = 1
        fd['stress'][2,0] = float(f[1])
        fd['stress'][2,1] = float(f[2])
        fd['stress'][2,2] = float(f[3])
      if(len(f) > 2 and f[0].upper() == "#M"):
        fd['mtemp'][f[1]] = f[2]      
      if(len(f) > 3 and f[0].upper() == "#C"):
        fd['c'][0] = int(f[1])
        fd['c'][1] = int(f[2])
        fd['c'][2] = int(f[3])
      if(len(f) > 1 and f[0].upper() == "#E"):
        fd['energy'] = float(f[1])
        fd['e'] = 1    
      if(len(f) > 1 and f[0].upper() == "#EPA"):
        fd['energy_per_atom'] = float(f[1])
        fd['e'] = 1   
      if(len(f) > 1 and f[0].upper() == "#RCUT"):
        fd['rcut'] = f[1]
        if(fd['rverlet'] == 0):
          fd['rverlet'] = f[1]        
      if(len(f) > 1 and f[0].upper() == "#RVERLET"):
        fd['rverlet'] = f[1]      
      if(len(f) > 1 and f[0].upper() == "#L_UNITS"):
        fd['l_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#E_UNITS"):
        fd['e_units'] = f[1]  
      if(len(f) > 1 and f[0].upper() == "#F_UNITS"):
        fd['f_units'] = f[1]
      if(len(f) > 1 and f[0].upper() == "#S_UNITS"):
        fd['s_units'] = f[1]
      if(len(f) >= 4 and f[0][0] != "#"):    
        count_coord = False
#print(f)
        if(len(f) >= 4):
          count_coord = True
        if(len(f) >= 7):
          fd['f'] = 1
        if(count_coord):
          fd['coord_count_prim'] = fd['coord_count_prim'] + 1

# COORD SIZE
    c = numpy.identity(3)
    
# EXPAND
    c[0,0] = float(fd['c'][0]) 
    c[1,1] = float(fd['c'][1]) 
    c[2,2] = float(fd['c'][2]) 
    fd['coord_count'] = (fd['c'][0] * fd['c'][1] * fd['c'][2]) * fd['coord_count_prim']
    fd['uv'] = numpy.matmul(c, fd['uv_prim'])
    
# MAKE ARRAYS
    fd['coords_label_prim'] = [None] * fd['coord_count_prim']
    fd['coords_label_id_prim'] = numpy.zeros((fd['coord_count_prim'],), dtype=numpy.int32,)
    fd['coords_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    if(fd['f'] == 1):
      fd['forces_prim'] = numpy.zeros((fd['coord_count_prim'],3,),)
    
    fd['coords_label'] = [None] * fd['coord_count']
    fd['coords_label_id'] = numpy.zeros((fd['coord_count'],), dtype=numpy.int32,)
    fd['coords'] = numpy.zeros((fd['coord_count'],3,),)
    if(fd['f'] == 1):
      fd['forces'] = numpy.zeros((fd['coord_count'],3,),)
        
    fd['n_atoms_prim'] = fd['coord_count_prim']
    fd['n_atoms'] = fd['coord_count']   
    
# READ COORDS
    n = 0
    for line in content:
      line = line.strip() 
      f = std.to_fields(line, ' ')
      if(len(f) >= 4 and f[0][0] != "#"):   
        if(len(f) >= 4):             
          label_str, label_id = labels.add(f[0])            
#print(f[0], label_str, label_id)
          fd['coords_label_prim'][n] = label_str
          fd['coords_label_id_prim'][n] = label_id
          fd['coords_prim'][n,0] = float(f[1])
          fd['coords_prim'][n,1] = float(f[2])
          fd['coords_prim'][n,2] = float(f[3])
          
        if(len(f) >= 7):
          fd['forces_prim'][n,0] = float(f[4])
          fd['forces_prim'][n,1] = float(f[5])
          fd['forces_prim'][n,2] = float(f[6])
        n = n + 1
    
# EXPAND COORDS
    m = 0
    for i in range(fd['c'][0]):
      for j in range(fd['c'][1]):
        for k in range(fd['c'][2]):
          for n in range(fd['coord_count_prim']):
            fd['coords_label'][m] = fd['coords_label_prim'][n]
            fd['coords_label_id'][m] = fd['coords_label_id_prim'][n]
            fd['coords'][m,0] = (i + fd['coords_prim'][n,0]) / fd['c'][0]
            fd['coords'][m,1] = (j + fd['coords_prim'][n,1]) / fd['c'][1]
            fd['coords'][m,2] = (k + fd['coords_prim'][n,2]) / fd['c'][2]
            if(fd['f'] == 1):
              fd['forces'][m,0] = fd['forces_prim'][n,0]
              fd['forces'][m,1] = fd['forces_prim'][n,1]
              fd['forces'][m,2] = fd['forces_prim'][n,2]
            m = m + 1

# SORT ENERGY
    if(fd['e'] == 1 and fd['energy_per_atom'] != None and fd['energy'] == None):
      try:
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass            
    if(fd['e'] == 1 and fd['energy'] != None and fd['energy_per_atom'] == None):
      try:
        fd['energy_per_atom'] = fd['energy'] / fd['n_atoms_prim']
        fd['energy'] = fd['energy_per_atom'] * fd['n_atoms']
      except:
        pass

# APPEND
    g.configs['configs'].append(fd)

#
#  FILE TYPE (standard, quantum espresso)
#
        
  @staticmethod
  def file_type(file_path):
    content = std.file_to_list(file_path)
    
# Check if standard file
    count = 0
    for line in content:
      fields = line.split(" ")
      if(fields[0].upper() == "#ALAT"):
        count = count + 1
      if(fields[0].upper() == "#X"):
        count = count + 1
      if(fields[0].upper() == "#Y"):
        count = count + 1
      if(fields[0].upper() == "#Z"):
        count = count + 1
    if(count >= 4):
      return 'std'
  
# Check if pwscf/qe file
    count = 0
    for line in content:
      if(line.strip()[0:13] == "Program PWSCF"):
        count = count + 1
      if(line.strip()[0:27] == "bravais-lattice index     ="):
        count = count + 1
      if(line.strip()[0:27] == "kinetic-energy cutoff     ="):
        count = count + 1
      if(line.strip()[0:27] == "mixing beta               ="):
        count = count + 1
      if(line.strip()[0:27] == "Exchange-correlation      ="):
        count = count + 1
      if(line.strip()[0:9] == "JOB DONE."):
        count = count + 1
    if(count >= 4):
      return 'qe'   
      
# Check if raw from pwscf/qe file
    count = 0
    for line in content:
      line = line.strip().upper()
      if(line == "RAW"):
        return 'raw'  

# Espresso files

  @staticmethod
  def qe(file_path):       
    qe = pwscf_output(file_path)
    xyz = qe.make_xyz()    
    n = 0
    for xyz_inner in xyz:
      n = n + 1
      configs.add_config(xyz_inner, file_path, n, config_type='QE')
      
  @staticmethod
  def raw(file_path):  

    fd = []    
    fh = open(file_path, 'r')
    for line in fh:
      line = std.one_space(line.strip())
      if(line != ""):
        fd.append(line)
    fh.close()
    
# Get atom count
    na = 0  
    i = 0
    while(i<len(fd)):
      if(fd[i][0:6].upper() == "#ATOMS"):
        na = int(fd[i+1])
        i = len(fd)
      i = i + 1
      
# Read data into vars/lists
    alat = ''
    energy = ''
    uv = []
    labels = []
    coords = []
    forces = []
    stress = []
    
    i = 0
    while(i<len(fd)):
      if(fd[i][0:5].upper() == "#ALAT"):
        i = i + 1
        alat = fd[i]
      elif(fd[i][0:3].upper() == "#UV"):
        for j in range(3):
          i = i + 1
          f = fd[i].split("(")
          f = f[2].split(")")
          f = f[0].strip().split(" ")
          uv.append(f)
      elif(fd[i][0:7].upper() == "#ENERGY"):
        i = i + 1
        f = fd[i].split("=")
        f = f[1].strip().split(" ")
        energy = f[0]
      elif(fd[i][0:7].upper() == "#COORDS"):
        for j in range(na):
          i = i + 1
          f = fd[i].strip().split(" ")
          labels.append(f[1])
          coords.append([f[6],f[7],f[8]])
      elif(fd[i][0:7].upper() == "#FORCES"):
        for j in range(na):
          i = i + 1
          f = fd[i].strip().split(" ")
          forces.append([f[6],f[7],f[8]])
      elif(fd[i][0:7].upper() == "#STRESS"):
        for j in range(3):
          i = i + 1
          f = fd[i].strip().split(" ")
          stress.append([f[3],f[4],f[5]])
          
# Increment
      i = i + 1
      
#print(na)
#print(alat)
#print(energy)
#print(uv)
#print(coords)
#print(forces)
    
    content = []
    content.append("#L_UNITS bohr")
    content.append("#E_UNITS ry")
    content.append("#F_UNITS ry/bohr")
    content.append("#S_UNITS kbar")
    content.append("#ALAT " + alat)
    content.append("#X " + uv[0][0] + " " + uv[0][1] + " " + uv[0][2])
    content.append("#Y " + uv[1][0] + " " + uv[1][1] + " " + uv[1][2])
    content.append("#Z " + uv[2][0] + " " + uv[2][1] + " " + uv[2][2])
    content.append("#SX " + stress[0][0] + " " + stress[0][1] + " " + stress[0][2])
    content.append("#SY " + stress[1][0] + " " + stress[1][1] + " " + stress[1][2])
    content.append("#SZ " + stress[2][0] + " " + stress[2][1] + " " + stress[2][2])
    content.append("#C 1 1 1")
    content.append("#E " + energy)
    content.append("#RCUT 6.5")
    
    for i in range(na):
      content.append(labels[i] + " " + coords[i][0] + " " + coords[i][1] + " " + coords[i][2] + " " + forces[i][0] + " " + forces[i][1] + " " + forces[i][2])
      
    configs.add_config(content, file_path, 1, config_type='RAW')

    print(content)

  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'configs.dat', 'w')
      
      n = 0
      for c in g.configs['configs']:
              
        n = n + 1
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('CONFIG ' + str(n) + '\n')
        fh.write('#############################################################\n')
        fh.write('#############################################################\n')
        fh.write('file_type  ' + str(c['file_type']) + ' \n')
        fh.write('file_path  ' + str(c['file_path']) + ' \n')
        fh.write('file_part  ' + str(c['file_part']) + ' \n')
        fh.write('alat  ' + str(c['alat']) + ' \n')
        fh.write('c         ' + str(c['c'][0])  + ' ' + str(c['c'][1]) + ' ' + str(c['c'][2]) + ' ' + ' \n')
        fh.write('uv_prim   ' + str(c['uv_prim'][0,0]) + ' ' + str(c['uv_prim'][0,1]) + ' ' + str(c['uv_prim'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][1,0]) + ' ' + str(c['uv_prim'][1,1]) + ' ' + str(c['uv_prim'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][2,0]) + ' ' + str(c['uv_prim'][2,1]) + ' ' + str(c['uv_prim'][2,2]) + ' ' + ' \n')
        fh.write('uv        ' + str(c['uv'][0,0]) + ' ' + str(c['uv'][0,1]) + ' ' + str(c['uv'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][1,0]) + ' ' + str(c['uv'][1,1]) + ' ' + str(c['uv'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][2,0]) + ' ' + str(c['uv'][2,1]) + ' ' + str(c['uv'][2,2]) + ' ' + ' \n')
        fh.write('stress    ' + str(c['stress'][0,0]) + ' ' + str(c['stress'][0,1]) + ' ' + str(c['stress'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][1,0]) + ' ' + str(c['stress'][1,1]) + ' ' + str(c['stress'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][2,0]) + ' ' + str(c['stress'][2,1]) + ' ' + str(c['stress'][2,2]) + ' ' + ' \n')
        fh.write('energy    ' + str(c['energy']) + ' \n')
        fh.write('epa       ' + str(c['energy_per_atom']) + ' \n')
        fh.write('rcut      ' + str(c['rcut']) + ' \n')
        fh.write('rverlet   ' + str(c['rverlet']) + ' \n')
        fh.write('e         ' + str(c['e']) + ' \n')
        fh.write('f         ' + str(c['f']) + ' \n')
        fh.write('s         ' + str(c['s']) + ' \n')    
        
        fh.write('\n')        
        fh.write('COORDS PRIM (' + str(c['coord_count_prim']) + ')\n')
        for k in range(c['coord_count_prim']):        
          fh.write(str(c['coords_label_prim'][k]) + '  ')  
          fh.write(str(c['coords_prim'][k,0]) + '  ')  
          fh.write(str(c['coords_prim'][k,1]) + '  ')  
          fh.write(str(c['coords_prim'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces_prim'][k,0]) + '  ')  
            fh.write(str(c['forces_prim'][k,1]) + '  ')  
            fh.write(str(c['forces_prim'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
        
        fh.write('\n')        
        fh.write('COORDS (' + str(c['coord_count']) + ')\n')
        for k in range(c['coord_count']):        
          fh.write(str(c['coords_label'][k]) + '  ')  
          fh.write(str(c['coords'][k,0]) + '  ')  
          fh.write(str(c['coords'][k,1]) + '  ')  
          fh.write(str(c['coords'][k,2]) + '  ')             
          if(c['f'] == 1):     
            fh.write(str(c['forces'][k,0]) + '  ')  
            fh.write(str(c['forces'][k,1]) + '  ')  
            fh.write(str(c['forces'][k,2]) + '  ')           
          fh.write('\n')  
        fh.write('\n')    
             
        fh.write('\n')               
        fh.write('\n')  
         
    fh.close()

  @staticmethod
  def save(cn, dir='', file_name=''):
  
    if(dir == ''):
      dir = g.dirs['configs']
    if(file_name == ''):
      file_name = str(cn + 1)
      while(len(file_name) < 4):
        file_name = "0" + file_name
      file_name = "config_" + file_name + ".dat"
  
    out = ""
    out += "/* CONFIG FILE */\n"
    out += "#L_UNITS ang\n"
    out += "#E_UNITS eV\n"
    out += "#S_UNITS GPA\n"
    out += "#F_UNITS EV/ANG\n"
    out += "#ALAT " + str(g.configs['configs'][cn]['alat']) + "\n"
    out += "#X " + str(g.configs['configs'][cn]['uv'][0,0]) + " "
    out += str(g.configs['configs'][cn]['uv'][0,1]) + " "
    out += str(g.configs['configs'][cn]['uv'][0,2]) + "\n"
    out += "#Y " + str(g.configs['configs'][cn]['uv'][1,0]) + " "
    out += str(g.configs['configs'][cn]['uv'][1,1]) + " "
    out += str(g.configs['configs'][cn]['uv'][1,2]) + "\n"
    out += "#Z " + str(g.configs['configs'][cn]['uv'][2,0]) + " "
    out += str(g.configs['configs'][cn]['uv'][2,1]) + " "
    out += str(g.configs['configs'][cn]['uv'][2,2]) + "\n"
    out += "#C 1 1 1\n"
    if(g.configs['configs'][cn]['s'] == 1):
      out += "#SX " + str(g.configs['configs'][cn]['uv'][0,0]) + " "
      out += str(g.configs['configs'][cn]['uv'][0,1]) + " "
      out += str(g.configs['configs'][cn]['uv'][0,2]) + "\n"
      out += "#SY " + str(g.configs['configs'][cn]['uv'][1,0]) + " "
      out += str(g.configs['configs'][cn]['uv'][1,1]) + " "
      out += str(g.configs['configs'][cn]['uv'][1,2]) + "\n"
      out += "#SZ " + str(g.configs['configs'][cn]['uv'][2,0]) + " "
      out += str(g.configs['configs'][cn]['uv'][2,1]) + " "
      out += str(g.configs['configs'][cn]['uv'][2,2]) + "\n"
    out += "#RCUT " + str(g.configs['configs'][cn]['rcut']) + "\n"
    if(g.configs['configs'][cn]['e'] == 1):
      out += "#EPA " + str(g.configs['configs'][cn]['energy_per_atom']) + "\n"
    for n in range(g.configs['configs'][cn]['coord_count']):
      out += str(g.configs['configs'][cn]['coords_label'][n]) + " "
      out += str(g.configs['configs'][cn]['coords'][n,0]) + " "
      out += str(g.configs['configs'][cn]['coords'][n,1]) + " "
      out += str(g.configs['configs'][cn]['coords'][n,2])
      if(g.configs['configs'][cn]['f'] == 1):      
        out += " "
        out += str(g.configs['configs'][cn]['forces'][n,0]) + " "
        out += str(g.configs['configs'][cn]['forces'][n,1]) + " "
        out += str(g.configs['configs'][cn]['forces'][n,2])
      out += "\n"
  
    fh = open(dir + "/" + file_name, 'w')
    fh.write(out)
    fh.close()
    
    """
    'file_type': '',  
    'file_path': '',  
    'file_part': 0,  
    'alat': 0.0,  
    'uv_prim': numpy.zeros((3,3,),),
    'uv': numpy.zeros((3,3,),),
    'stress': numpy.zeros((3,3,),),
    'c': numpy.zeros((3,),dtype=numpy.int32,),
    'rcut': 0.0,
    'rverlet': 0.0,
    'mtemp': {},
    'm': {},
    'coord_count_prim': 0,
    'coords_label_prim': None,
    'coords_label_id_prim': None,
    'coords_prim': None,
    'forces_prim': None,
    'coord_count': 0,
    'coords_label': None,
    'coords_label_id': None,
    'coords': None,
    'forces': None,
    'energy_per_atom': None,
    'energy': None,
    'e': 0,
    'f': 0,
    's': 0,
    'n_atoms_prim': 0,
    'n_atoms': 0,
    'l_units': 'ANG',
    'e_units': 'EV',
    'f_units': 'EV/ANG',
    's_units': 'EV/ANG3',
    """

# Run after configs read in
# After dft_energy_adjustments
    
  @staticmethod
  def complete(): 
#print("Complete Configs")
#print(g.configs['configs'])
    
# CONVERT
    for i in range(len(g.configs['configs'])):
    
      g.configs['configs'][i]['alat'] = units.convert(g.configs['configs'][i]['l_units'], 'ang', g.configs['configs'][i]['alat'])
      g.configs['configs'][i]['l_units'] = 'ANG'
      
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
      g.configs['configs'][i]['energy_per_atom'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy_per_atom'])
      g.configs['configs'][i]['energy'] = units.convert(g.configs['configs'][i]['e_units'], 'ev', g.configs['configs'][i]['energy'])
      g.configs['configs'][i]['e_units'] = 'EV'
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
      
# print(g.configs['configs'][i]['f_units'])
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces_prim'])):
          for nb in range(len(g.configs['configs'][i]['forces_prim'][na])):
            g.configs['configs'][i]['forces_prim'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces_prim'][na,nb])
      
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['forces'])):
          for nb in range(len(g.configs['configs'][i]['forces'][na])):
            g.configs['configs'][i]['forces'][na,nb] = units.convert(g.configs['configs'][i]['f_units'], 'EV/ANG', g.configs['configs'][i]['forces'][na,nb])
     
      if(g.configs['configs'][i]['f'] == 1):
        for na in range(len(g.configs['configs'][i]['stress'])):
          for nb in range(len(g.configs['configs'][i]['stress'][na])):
            g.configs['configs'][i]['stress'][na,nb] = units.convert(g.configs['configs'][i]['s_units'], 'EV/ANG3', g.configs['configs'][i]['stress'][na,nb])
    
# ADJUST ENERGY QE
    for i in range(len(g.configs['configs'])):
      if(g.configs['configs'][i]['file_type'] == 'QE'):
#print(g.configs['configs'][i]['coord_count'], g.configs['configs'][i]['coord_count_prim'])
        for k in range(g.configs['configs'][i]['coord_count']):  
          e_adj = g.dft_energy_adjustments[g.configs['configs'][i]['coords_label_id'][k]]          
#except:
#  print("Energy Adjustment Missing for " + str(k))
#  eampa.exit()
          g.configs['configs'][i]['energy'] = g.configs['configs'][i]['energy'] + e_adj['calc_apaev']
        g.configs['configs'][i]['energy_per_atom'] = g.configs['configs'][i]['energy'] / g.configs['configs'][i]['coord_count'] 
#print(g.configs['configs'][i]['energy'], g.configs['configs'][i]['energy_per_atom'])
          
#print(e_adj['calc_apaev'])

#c['coords_label'][k]
#print('QE')
        
#  g.dft_energy_adjustments[label_id]
  
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def efs_add_config():  

    while(len(g.configs['config_results'])<len(g.configs['configs'])):
      g.configs['config_results'].append({'f_id': None, 'energy': None, 'stress': None, 'force': None,})

    for n in range(len(g.configs['configs'])):
      c = g.configs['configs'][n]

#print(c)
      efs.add_config( 
                     6.5,
                     c['alat'], 
                     c['uv'], 
                     c['coords_label_id'], 
                     c['coords'], 
                     c['energy'], 
                     c['forces'], 
                     c['stress'], 
                     c['f'], 
                     c['s']
                    )    
      g.configs['configs'][n]['f_id'] = efs.cc
      g.configs['config_results'][n]['f_id'] = efs.cc
    efs.make_nl()

  @staticmethod
  def efs_results(): 
    while(len(g.configs['config_results'])<len(g.configs['configs'])):
      cc = g.configs['config_results'][n]['f_id']
      print(efs.energies[cc])

#@staticmethod
#def efs_add_config():
    
###########################################
#  CLASS pwscf_outpu
###########################################
class pwscf_output:

  def __init__(self, file_in=None):
    self.reset()
    if(file_in != None):
      self.load(file_in)
      self.calcs()

  def reset(self):
  
    self.z = numpy.zeros((3,3))
    
# Important so store in it's own variable
    self.atom_count = 1   

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "error_code": None,
      "converged": False,
      "converged_in": None,
      
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "species": {},
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
      
      "mass_per_crystal": 0.0,
      "density": [],
    }
    
# Defaults
    self.xyz_units = 'evang'
    self.stress_units = 'gpa'
    
# Constants
    self.avogadro = 6.02214086E23

  def scf_settings(self):
    return {
    "bravais_lattice_index": None,
    "alat": None,
    "volume": None,
    "electrons": None, 
    "electrons_up": None, 
    "electrons_down": None, 
    "ecut_wfc": None, 
    "ecut_rho": None, 
    "convergence_threshold": None, 
    "mixing_beta": None, 
    "atomic_species": {},
    }
    
  def scf_crystal(self):
    return {
    "alat": 0.0,
    "cell_parameters": numpy.zeros((3,3)),
    "position_units": None,
    "atomic_labels": [],
    "atomic_positions": numpy.zeros((self.atom_count,3)),    
    "alat_adj": 0.0,
    "cell_parameters_adj": numpy.zeros((3,3)),
    "crystal_positions": numpy.zeros((self.atom_count,3)),
    "cell_volume": 0.0,
    "cell_density": 0.0,
    }

  def scf_results(self):
    return {
    "energy": 0.0,
    "total_force": 0.0,
    "stress": numpy.zeros((3,3)),
    "forces": numpy.zeros((self.atom_count,3)),
    "f_on": False,
    "s_on": False,
    }

#  Load, and use in another program
  def load(self, file_name): 
  
# Load data from file
    data_file = self.load_from_file(file_name)
    self.d = data_file.split("\n") 
  
#print(self.d)
  
# Reset data store
    self.reset()
       
# Load
    self.load_status()
    self.load_type()
    self.load_times()
    self.load_count()
    self.load_cpuinfo()
    self.load_crystal()
    self.load_scf_settings()
    self.load_results()
    self.load_species()
    
#print(len(self.data['crystals']))
#print(len(self.data['results']))
#self.load_scf('final_scf')
#self.load_results('initial_scf')
#self.load_results('final_scf')
    
  def load_status(self):    
# OK
###################################
    self.data['ok'] = False
    counter = 0

    for line in self.d:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        self.data['job_done'] = True
      if(pwscf_output.compare(line, "Exit code:")):
        self.data['error'] = True              
        try:
          self.data['error_code'] = int(line[13:-1])
        except:
          pass
      if(pwscf_output.compare(line, "convergence has been achieved")):
        self.data['converged'] = True
        try:
          self.data['converged_in'] = int(line[39:42])
        except:
          pass
    if(self.data['job_done']):
      self.data['ok'] = True
  
  def load_type(self):
# Calc Type
###################################
    self.data['type'] = "SCF"
    for line in self.d:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
    
  def load_times(self):
# Calc Type
###################################
    for line in self.d:
      line = line.strip()
      if(line[0:14] == "PWSCF        :"):
        fa = line.split(':')
        fb = fa[1].split('CPU')
        cpu = fb[0].strip()
        fc = fb[1].split('WALL')
        wall = fc[0].strip()   
        self.data['cpu_time'] = pwscf_output.time_in_seconds(cpu)  
        self.data['wall_time'] = pwscf_output.time_in_seconds(wall)  
#print(self.data['cpu_time'])
#print(self.data['wall_time'])
      
  @staticmethod
  def time_in_seconds(inp): 
    h = 0.0
    m = 0.0
    s = 0.0
    if('h' in inp):
      f = inp.split('h')
      try:  
        h = float(f[0])
      except:
        pass
      inp = f[1]
    if('m' in inp):
      f = inp.split('m')
      try:  
        m = float(f[0])
      except:
        pass
      inp = f[1]
    if('s' in inp):
      f = inp.split('s')
      try:  
        s = float(f[0])
      except:
        pass
    return 3600 * h + 60 * m + s
        
    """    
PWSCF        : 15m15.00s CPU    16m33.12s WALL 
PWSCF        :    28.10s CPU        29.31s WALL 
PWSCF        :     1h29m CPU        1h32m WALL         
    """
   
  def load_count(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "number of atoms/cell      =")):
        count = pwscf_output.extract(line, "=", "", "i")  
        try:
          self.atom_count = int(count)
          return self.atom_count 
        except:
          return 0
          
  def load_cpuinfo(self):
    counter = 0
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
          elif(pwscf_output.compare(line, "Threads/MPI process:")):
            self.data['threads_per_mpi_process'] = pwscf_output.extract(line, ":", "", "i") 
   
  def load_species(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line[0:14] == "atomic species"):
        n, line, line_uc = self.next_line(n, self.d)
        while(line.strip() != ""):      
          line = pwscf_output.single_spaces(line).strip()
          f = line.split(" ")
          self.data['scf_settings'][f[0]] = [float(f[1]), float(f[2])]   # Valence electrons, atomic mass
          n, line, line_uc = self.next_line(n, self.d)
        break
      
#self.data['scf_settings']
#
#  atomic species   valence    mass     pseudopotential
#  Al             3.00    26.98200     Al( 1.00)
 
  def load_scf_settings(self):
    self.data['scf_settings'] = self.scf_settings()  
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "bravais-lattice index     =")):
        self.data['scf_settings']['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
      elif(pwscf_output.compare(line, "lattice parameter (alat)  =")):
        self.data['scf_settings']['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
      elif(pwscf_output.compare(line, "unit-cell volume          =")):
        self.data['scf_settings']['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
      elif(pwscf_output.compare(line, "number of atoms/cell      =")):
        self.data['scf_settings']['nat'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of atomic types    =")):
        self.data['scf_settings']['types'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of electrons       =")):
        str_e = pwscf_output.extract(line, "=", "", "s")
        e, eu, ed = pwscf_output.electron_string(str_e)
        self.data['scf_settings']['electrons'] = e
        self.data['scf_settings']['electrons_up'] = eu
        self.data['scf_settings']['electrons_down'] = ed
      elif(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
        self.data['scf_settings']['ecut_wfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "charge density cutoff     =")):
        self.data['scf_settings']['ecut_rho'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "convergence threshold     =")):
        self.data['scf_settings']['convergence_threshold'] = pwscf_output.extract(line, "=", "", "f")  
      elif(pwscf_output.compare(line, "mixing beta               =")):
        self.data['scf_settings']['mixing_beta'] = pwscf_output.extract(line, "=", "", "f")  
      elif(("atomic species" in line) and ("valence" in line) and ("mass" in line) and ("pseudopotential" in line)):
        loop = True
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)
          if(line.strip() == ""):
            loop = False
          else:  
            line = pwscf_output.single_spaces(line).strip()
            line_arr = line.split(" ")
            self.data['scf_settings']['atomic_species'][line_arr[0]] = {}
            self.data['scf_settings']['atomic_species'][line_arr[0]]['valence'] = line_arr[1]
            self.data['scf_settings']['atomic_species'][line_arr[0]]['mass'] = line_arr[2]

# End of file/loop
        n = len(self.d)

###################################
# LOAD CRYSTALS FROM OUTPUT FILE
###################################

  def load_crystal(self):
# Make new list for crystals
    self.data['crystals'] = []
    
# FIRST
    n = 0
    crystal = self.scf_crystal()
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)  
      if(line[0:10] == "celldm(1)="):
        crystal['alat'] = float(line[10:21].strip())
      elif(pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(pwscf_output.compare(line.strip(), "Cartesian axes")):
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
# Unit
        crystal['position_units'] = "alat"
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
        
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
        
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# MIDDLE
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:15] == "CELL_PARAMETERS"):
# Create
        crystal = self.scf_crystal()
        
# Get alat
        line_arr = line.split("=")
        line_arr = line_arr[1].split(")")
        crystal['alat'] = float(line_arr[0].strip())
        
#Cell Parameters
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          line = pwscf_output.single_spaces(line)
          fields = line.split(" ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(line[0:9] == "density ="):  
        crystal['density'] = float(line[10:23])      
            
      elif(line[0:16] == "ATOMIC_POSITIONS"):     
        
# Unit
        crystal['position_units'] = "crystal"
        
# Read Coords
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          elif(line.strip() == "End final coordinates"):
            loop = False
          else:
            line = pwscf_output.single_spaces(line)
            line_arr = line.split(" ")
            crystal['atomic_labels'].append(line_arr[0])
            
            crystal['atomic_positions'][k, 0] = float(line_arr[1])
            crystal['atomic_positions'][k, 1] = float(line_arr[2])
            crystal['atomic_positions'][k, 2] = float(line_arr[3]) 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
    
# END
    n = 0
    crystal = self.scf_crystal()
    d = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:10] == "celldm(1)="):
        d = d + 1
        if(d == 2):
          crystal['alat'] = float(line[10:21].strip())
      elif(d == 2 and pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
     
      elif(d == 2 and pwscf_output.compare(line.strip(), "Cartesian axes")):      
        
# Unit
        crystal['position_units'] = "alat"
        
# Read coords
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
# Increment
            k = k + 1
            
# Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
# Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
# Loop through crystals
    for i in range(len(self.data['crystals'])):
# Adjust alat and cell_parameters so celldm(1)=1.0
      factor = 1.0 / self.data['crystals'][i]['cell_parameters'][0, 0]
      
      self.data['crystals'][i]['alat_adj'] = self.data['crystals'][i]['alat'] * self.data['crystals'][i]['cell_parameters'][0, 0]
      self.data['crystals'][i]['cell_parameters_adj'][:, :] = factor * self.data['crystals'][i]['cell_parameters'][:, :]
    
# Make crystal_positions
      if(self.data['crystals'][i]['position_units'] == 'crystal'):
        self.data['crystals'][i]['crystal_positions'][:,:] = self.data['crystals'][i]['atomic_positions'][:,:] 
      elif(self.data['crystals'][i]['position_units'] == 'alat'):  
        minv = numpy.linalg.inv(self.data['crystals'][i]['cell_parameters'][:, :])
        for j in range(len(self.data['crystals'][i]['atomic_positions'])):
          self.data['crystals'][i]['crystal_positions'][j, :] = numpy.matmul(minv[:,:], self.data['crystals'][i]['atomic_positions'][j, :])

  def load_results(self):
  
# Make new list for results
    self.data['results'] = []
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      
# READ ENERGY
      if(pwscf_output.compare(line, "!    total energy")):
# Create dictionary
        results = self.scf_results()  
        results['energy'] = pwscf_output.extract(line, "=", "Ry", "f") 
        
# READ FORCES
      elif(pwscf_output.compare(line, "Forces acting on atoms")):
        n, line, line_uc = self.next_line(n, self.d)  
        loop = True
        f = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)  
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("force =")
            fields = pwscf_output.single_spaces(line_arr[1].strip()).split(" ")
            results['forces'][f,0] = float(fields[0])
            results['forces'][f,1] = float(fields[1])
            results['forces'][f,2] = float(fields[2])
            f = f + 1
        if(f>0):
          results['f_on'] = True
        
# READ TOTAL FORCE
      elif(pwscf_output.compare(line, "Total force =")):
        results['total_force'] = pwscf_output.extract(line, "=", "T", "f")
        
# READ STRESS
      elif(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):  
        results['s_on'] = True      
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)  
          fields = pwscf_output.extract(line, "", "", "f", " ", True)  
          results['stress'][j,0] = fields[3] 
          results['stress'][j,1] = fields[4] 
          results['stress'][j,2] = fields[5]

#SAVE
        self.data['results'].append(results)

  def calcs(self):
  
# MASS PER CELL
    self.data['mass_per_crystal'] = 0.0    
    for l in self.data['crystals'][0]['atomic_labels']:
      self.data['mass_per_crystal'] = self.data['mass_per_crystal'] + self.data['scf_settings'][l][1]    
      
# CELL VOLUMES
    for i in range(len(self.data['crystals'])):
      self.data['crystals'][i]['cell_volume'] = pwscf_output.cell_volume(self.data['crystals'][i]['alat'], self.data['crystals'][i]['cell_parameters'][:,:])  
      self.data['crystals'][i]['cell_density'] = pwscf_output.cell_density(self.data['crystals'][i]['cell_volume'], self.data['mass_per_crystal'])
      
  def next_line(self, n, data):
    if(n < len(data)):
      line = data[n].strip()
      line_uc = line.upper()
      n = n + 1
      return n, line, line_uc
    else:
      n = n + 1
      return n, None, None
    
  def store(self, store, line, field, n=0):
    l, f = pwscf_output.read_line(line, field)  
    if(l != False):
      self.data[store] = f[n]

#  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0

#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
# Init variable
    file_data = ""

# Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data

#################################
# Get
#################################

  def get_nat(self): 
    return self.atom_count

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['scf_settings']['volume'] 
    
  def get_volume_per_atom(self):
    return self.data['scf_settings']['volume'] / self.atom_count
  
  def get_total_energy(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy']  
    else:
      return self.data['results'][n]['energy']  
    
  def get_energy_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy'] / self.atom_count  
    else:
      return self.data['results'][n]['energy'] / self.atom_count 
    
  def get_total_force(self, n = None):
    if(n == None):
      return self.data['results'][-1]['total_force']  
    else:
      return self.data['results'][n]['total_force']
    
  def get_force_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['total_force'] / self.atom_count  
    else:
      return self.data['results'][n]['total_force'] / self.atom_count  
  
  def get_density(self):
    return self.data['density']  
    
  def get_cell_parameters(self):
    cp = ['alat', 
          [str(self.data['crystal_calc'][0,0]), str(self.data['crystal_calc'][0,1]), str(self.data['crystal_calc'][0,2])], 
          [str(self.data['crystal_calc'][1,0]), str(self.data['crystal_calc'][1,1]), str(self.data['crystal_calc'][1,2])], 
          [str(self.data['crystal_calc'][2,0]), str(self.data['crystal_calc'][2,1]), str(self.data['crystal_calc'][2,2])]]
    return cp

# Return relaxed unit vector
  def get_cell_array(self):
    return self.data['crystal_calc']

# return alat and normalised unit vector
  def get_norm_relaxed(self):
    alat = self.data['crystals'][-2]['alat_adj']
    cp = self.data['crystals'][-2]['cell_parameters_adj']
    return alat, cp

# Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  
  def get_job_done(self):
    return self.data['job_done']
    
  def get_job_error(self):
    return self.data['error']
    
  def get_job_converged(self):
    return self.data['converged']

  def get_ok(self):
    return self.data['ok']

# GET - VC-RELAXED

  def get_vc_relax_alat(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['alat_adj']
    else:
      return None

  def get_vc_relax_cp(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_parameters_adj']
    else:
      return None

  def get_vc_relax_volume(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_volume']
    else:
      return None

  def get_vc_relax_density(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_density']
    else:
      return None
      
  def get_mass_per_crystal(self):
    return self.data['mass_per_crystal']
      
  def get_times(self):
    return self.data['cpu_time'], self.data['wall_time']
      
#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()

  def print_menu(self):
    pwscf_output.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")

  def i_load(self):
    pwscf_output.header("Load Output File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    print("File loaded.")
    input()

  def i_display(self):
    pwscf_output.header("Display File")
    self.output_details()
    input()

  def output_details(self):
    print("Output")
    print("=======================================================================")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    print("=======================================================================")
    print()
    
  def xyz_evang(self):
    self.xyz_units = 'ev/ang'
    self.stress_units = 'gpa'

  def xyz_stress_gpa(self):
    self.stress_units = 'gpa'

  def make_xyz(self, option=None):
    self.xyz = []
    if(option == None):
      for rn in range(len(self.data['results'])):
        option = rn + 1
    elif(option == -1):
      option = len(self.data['results'])
    else:
      option = (option - 1) % len(self.data['results']) + 1
    self.make_xyz_inner(option)
    return self.xyz
    
  def make_xyz_inner(self, option):
    if(len(self.data['results'])==0):
      return False
    if(len(self.data['crystals'])==0):
      return False  
   
    rn = (option - 1) % len(self.data['results'])
    cn = rn
    if(rn == len(self.data['results']) - 1):
      cn = len(self.data['crystals']) - 1
      
    crystal = self.data['crystals'][cn]
    result = self.data['results'][rn]
    settings = self.data['scf_settings']
    species = settings['atomic_species']
    
# Add new list and set counter n
    self.xyz.append([])
    n = len(self.xyz) - 1
    
# Add data
    self.xyz[n].append("#ALAT " + str(crystal['alat_adj']))
    self.xyz[n].append("#X " + str(crystal['cell_parameters_adj'][0][0]) + " " + str(crystal['cell_parameters_adj'][0][1]) + " " + str(crystal['cell_parameters_adj'][0][2]))
    self.xyz[n].append("#Y " + str(crystal['cell_parameters_adj'][1][0]) + " " + str(crystal['cell_parameters_adj'][1][1]) + " " + str(crystal['cell_parameters_adj'][1][2]))
    self.xyz[n].append("#Z " + str(crystal['cell_parameters_adj'][2][0]) + " " + str(crystal['cell_parameters_adj'][2][1]) + " " + str(crystal['cell_parameters_adj'][2][2]))
    
# Just use 1 1 1
    self.xyz[n].append("#C 2 2 2")
    self.xyz[n].append("#RCUT 6.5")
    self.xyz[n].append("#L_UNITS bohr")
    self.xyz[n].append("#E_UNITS ry")
    self.xyz[n].append("#F_UNITS ry/bohr")
    self.xyz[n].append("#S_UNITS kbar")
    
    self.xyz[n].append("#E " + str(result['energy']))
     
    if(result['s_on']):
      self.xyz[n].append("#SX " + str(result['stress'][0,0]) + " " +  str(result['stress'][0,1]) + " " + str(result['stress'][0,2]))
      self.xyz[n].append("#SY " + str(result['stress'][1,0]) + " " +  str(result['stress'][1,1]) + " " + str(result['stress'][1,2]))
      self.xyz[n].append("#SZ " + str(result['stress'][2,0]) + " " +  str(result['stress'][2,1]) + " " + str(result['stress'][2,2]))
      
    for label in species.keys():
      self.xyz[n].append("#M " + label + " " + str(species[label]['mass']))
    
    for i in range(self.atom_count):
      line = crystal['atomic_labels'][i]
      line = line + " " + str(crystal['crystal_positions'][i,0]) + " " + str(crystal['crystal_positions'][i,1]) + " " + str(crystal['crystal_positions'][i,2])      
      if(result['f_on']):   
        line = line + " " + str(result['forces'][i,0]) + " " + str(result['forces'][i,1]) + " " + str(result['forces'][i,2])
#"s_on": False,
      self.xyz[n].append(line)
      
# Return
    return self.xyz[n]
    
  def get_data(self, file=None): 
  
    """
# Important so store in it's own variable
    self.atom_count = 1    

# Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
    }
    """
    
    out = "##############################################################################################################\n"
    out = out + "atom count:                  " + str(self.atom_count) + "\n"
    out = out + "ok:                          " + str(self.data['ok']) + "\n"
    out = out + "job_done:                    " + str(self.data['job_done']) + "\n"
    out = out + "error:                       " + str(self.data['error']) + "\n"
    out = out + "type:                        " + str(self.data['type']) + "\n"
    out = out + "summary:                     " + str(self.data['summary']) + "\n"
    out = out + "mpi_processes:               " + str(self.data['mpi_processes']) + "\n"
    out = out + "threads_per_mpi_process:     " + str(self.data['threads_per_mpi_process']) + "\n"
    out = out + "scf_settings:                " + str(len(self.data['scf_settings'])) + "\n"
    
    for k in self.data['scf_settings'].keys():    
      out = out + "                             " + k + '  '+ str(self.data['scf_settings'][k]) + "\n"  
    
    out = out + "crystals (count):            " + str(len(self.data['crystals'])) + "\n"
    out = out + "results (count):             " + str(len(self.data['results'])) + "\n"
    i = 0
    for i in range(len(self.data['crystals'])-1):
      out = out + "##############################################################################################################\n"
      c = self.data['crystals'][i]
      n = i + 1
      if(n == 1):
        out = out + "crystal " + str(n) + " (input):\n"         
      elif(n == len(self.data['crystals']) - 1):
        out = out + "crystal " + str(n) + " (relaxed):\n" 
      else:
        out = out + "crystal " + str(n) + ":\n"  
      out = out + "##############################################################################################################\n"
      out = out + "alat:                        " + str(c['alat']) + "\n"  
      out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
      out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
      out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 
      if(i<len(self.data['crystals'])-2):
        s = self.data['results'][i]
        out = out + "energy:                      " + str(s['energy']) + "\n"  
        out = out + "total_force:                 " + str(s['total_force']) + "\n"  
        out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
        out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
        for m in range(self.atom_count):
          out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
      
      out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
      out = out +     "density:                     " + str(c['cell_density']) + "\n" 
       
    c = self.data['crystals'][-1]
    s = self.data['results'][-1]
    n = len(self.data['crystals'])
    
    out = out + "##############################################################################################################\n"
    out = out + "crystal " + str(n) + " (final): \n"     
    out = out + "##############################################################################################################\n"
    out = out + "alat:                        " + str(c['alat']) + "\n" 
    out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
    out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
    out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 

    out = out + "energy:                      " + str(s['energy']) + "\n"  
    out = out + "total_force:                 " + str(s['total_force']) + "\n"  
    out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
    out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
    for m in range(self.atom_count):
      out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
         
    out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
    out = out +     "density:                     " + str(c['cell_density']) + "\n" 
    
    fh = open(file,'w')
    fh.write(out)
    fh.close()
    
    return out
    
#################################
# Static Methods
#################################

  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
    
  @staticmethod
  def extract(input_string, start=None, end=None, type=None, split=None, trim=False):
    if(start == ""):
      start = None
    if(end == ""):
      end = None
    if(trim):
      input_string = input_string.strip()
    
# Start/End
    start_n = None
    end_n = None
      
    if(start == None and end == None):   
      start_n = 0
      end_n = len(input_string)
    elif(start == None and end != None):  
      end_l = len(end)   
      start_n = 0
      for n in range(len(input_string)):
        if(input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
    elif(start != None and end == None):  
      start_l = len(start)
      end_n = len(input_string)
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
    else:  
      start_l = len(start)
      end_l = len(end)  
    
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
        if(start_n != None and input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
        
# Read
    result = input_string[start_n:end_n].strip()       

# Split
    if(split != None):
      if(split == " "):
#result = re.sub(r'\s\s+', ' ', result)
        result = pwscf_output.single_spaces(result)
      result = result.split(split)
      for i in range(len(result)):
        if(type.lower() == "f"):
          result[i] = float(result[i])
        elif(type.lower() == "i"):
          result[i] = int(result[i])
        
    else:  
      if(type.lower() == "f"):
        result = float(result)
      elif(type.lower() == "i"):
        result = int(result)
        
# Return
    return result
      
  @staticmethod
  def compare(line, field):
    line = line.strip()
    line = line.upper() 
    
    field = field.strip()
    field = field.upper()
    
    f_len = len(field)
    if(len(line) >= f_len and line[0:f_len] == field[0:f_len]):
      return True
    return False
    
  @staticmethod
  def read_line(line, field):
    line = line.strip()
#line = re.sub(r'\s\s+', ' ', line)
#line = re.sub(r'\s=\s', '=', line)
    line = pwscf_output.clean(line)
    line_uc = line.upper() 
    
    field = field.strip()
#field = re.sub(r'\s\s+', ' ', field)
#field = re.sub(r'\s=\s', '=', field)
    field = pwscf_output.clean(field)
    field = field.upper()
    
    f_len = len(field)
    if(len(line_uc) >= f_len and line_uc[0:f_len] == field[0:f_len]):
      output = line[f_len:].strip()
      fields = output.split(" ")
      return output, fields      
    return False, False
    
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False

  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

  @staticmethod
  def header(sub_title=""):
    pwscf_output.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_output.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  

  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
  @staticmethod
  def is_zero(arr):
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]):
        if(arr[i, j] != 0.0):
          return False
    return True
    
  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
# Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
    
# Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
# Add to string
      if(ok):
        str_out += char
    return str_out    
    
  @staticmethod
  def electron_string(str_in):
    arr = str_in.split("(up:")
    e = arr[0]
    if(len(arr) == 1):
      return e.strip(), None, None
    if(len(arr)==2):
      arr_b = arr[1].split(", down:")
      eu = arr_b[0]
      arr_c = arr_b[1].split(")")
      ed = arr_c[0]
      return e.strip(), eu.strip(), ed.strip()
  
    print("TEST")
    return "","",""
    
  @staticmethod
  def single_spaces(str_in):
    str_out = ""
    last = None
    for char in str_in:
      if(char != " " or (char == " " and last != " ")):
        str_out = str_out + char
      last = char
    return str_out
    
  @staticmethod
  def cell_volume(alat, cp):
    cp_alat = numpy.zeros((3,3,),)
    cp_alat[:,:] = alat * cp[:,:]
    v = numpy.dot(cp_alat[0,:],numpy.cross(cp_alat[1,:], cp_alat[2,:]))
    return v
    
  @staticmethod
  def cell_density(v, mpc):
    v_m3 = v * 1.48036e-31
    m = mpc * 1.66054E-027
    rho = m / v_m3
    return rho
    
"""
        
  def aaa():
    
# Load
###################################
    n = 0
    counter = 0
    while(n < len(data)):
      n, line, line_uc = self.next_line(n, data)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
        
          if(pwscf_output.compare(line, "bravais-lattice index     =")):
            self.data['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
            
          if(pwscf_output.compare(line, "lattice parameter (alat)  =")):
            self.data['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
            
          if(pwscf_output.compare(line, "unit-cell volume          =")):
            self.data['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
            
          if(pwscf_output.compare(line, "number of atoms/cell      =")):
            self.data['nat'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of atomic types    =")):
            self.data['types'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of electrons       =")):
            str_e = pwscf_output.extract(line, "=", "", "s")
            e, eu, ed = pwscf_output.electron_string(str_e)
            self.data['electrons'] = e
            self.data['electrons_up'] = eu
            self.data['electrons_down'] = ed
          
          if(pwscf_output.compare(line, "number of Kohn-Sham states=")):
            self.data['ks_states'] = pwscf_output.extract(line, "=", "", "i")   
            
          if(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
            self.data['ecutwfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
            
          if(pwscf_output.compare(line, "charge density cutoff     =")):
            self.data['ecutrho'] = pwscf_output.extract(line, "=", "Ry", "f")   
        
          if(pwscf_output.compare(line.strip(), "crystal axes:") and pwscf_output.is_zero(self.data['crystal_in'])):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_in'][j,:] = fields  
              self.data['crystal_calc'][j,:] = fields  
          
          if(pwscf_output.compare(line.strip(), "crystal axes:")):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_calc'][j,:] = fields            
        
          if(pwscf_output.compare(line, "!    total energy")):
            self.data['total_energy'] = pwscf_output.extract(line, "=", "Ry", "f")
            
          if(pwscf_output.compare(line, "Total force =")):
            self.data['total_force'] = pwscf_output.extract(line, "=", "T", "f")
            
          if(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):        
            self.data['stress_sum'] = 0.0
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)   
              fields = pwscf_output.extract(line, "", "", "f", " ", True)  
              self.data['stress'][j,0] = fields[3] 
              self.data['stress'][j,1] = fields[4] 
              self.data['stress'][j,2] = fields[5]
              self.data['stress_sum'] = self.data['stress_sum'] + abs(fields[0]) + abs(fields[1]) + abs(fields[2])
            
#                  "stress": numpy.zeros((3,3)),
#      "stress_sum": None,
            
          if(pwscf_output.compare(line, "density = ")):
            self.data['density_full'] = pwscf_output.extract(line, "=", "", "s")
            self.data['density'] = pwscf_output.extract(line, "=", "g/cm^3", "f")
          
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['cpu_time'] = pwscf_output.extract(line, ":", "CPU", "s")
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['wall_time'] = pwscf_output.extract(line, "CPU", "WALL", "s")
          
          if(pwscf_output.compare(line, "JOB DONE.")):
            self.data['job_done'] = True
            
          if(pwscf_output.compare(line, "Exit code:")):
            self.data['error'] = True  
"""

###########################################
#  CLASS b_prop
###########################################
class b_props:

  @staticmethod
  def load():
    if('bp' not in g.inp.keys()):
      return None
    if('bp_file' not in g.inp['bp'].keys()):
      return None
    try: 
      dir = g.inp['bp']['dir'].strip()
      bp_file = std.path(g.inp['bp']['dir'], g.inp['bp']['bp_file'])
    except:
      dir = ""
      bp_file = g.inp['bp']['bp_file']

    if(not(os.path.isfile(bp_file))):
      return None
    
# Read BP data file
    bp_inp = read_config.read_file(bp_file)

# READ IN UNITS
    try:
      bp_pressure = bp_inp['units']['pressure']
    except:
      bp_pressure = 'GPA'
    try:
      bp_length = bp_inp['units']['length']
    except:
      bp_length = 'ang'
    try:
      bp_energy = bp_inp['units']['energy']
    except:
      bp_energy = 'ev'
      
#print(bp_inp)
    
    for k in bp_inp.keys():
      if('potlabel' in bp_inp[k].keys() and 'alat' in bp_inp[k].keys()):
        potlabel = bp_inp[k]['potlabel']
        label_str, label_id = labels.add(potlabel)     
        
        newbp = b_props.make(label_id, label_str)
        newbp['bp_key_str'] = k

# There must be an alat value set
        newbp['alat'] = float(bp_inp[k]['alat'])
        newbp['alat'] = units.convert(bp_length, 'ang', newbp['alat'])
#newbp['label_str'] = label_str
                
        try:
          if(bp_inp[k]['type'].lower() == 'sc'):
            newbp['type'] = 1
            newbp['type_text'] = 'Simple Cubic'
          elif(bp_inp[k]['type'].lower() == 'bcc'):
            newbp['type'] = 2
            newbp['type_text'] = 'Body Centered Cubic'
          elif(bp_inp[k]['type'].lower() == 'fcc'):
            newbp['type'] = 3
            newbp['type_text'] = 'Face Centered Cubic'
          elif(bp_inp[k]['type'].lower() == 'zb'):
            newbp['type'] = 4
            newbp['type_text'] = 'Zinc Blende'
        except:
          pass 
          
        newbp['uv'][:,:] = 0.0
        newbp['uv'][0,0] = 1.0
        newbp['uv'][1,1] = 1.0
        newbp['uv'][2,2] = 1.0
        
        try: 
          if('uv' in bp_inp[k].keys()):
# CUBIC
            if(len(bp_inp[k]['uv']) == 1):
              newbp['uv'][0,0] = float(bp_inp[k]['uv'][0])
              newbp['uv'][1,1] = float(bp_inp[k]['uv'][0])
              newbp['uv'][2,2] = float(bp_inp[k]['uv'][0])
# CUBIC
            if(len(bp_inp[k]['uv']) == 3):
              newbp['uv'][0,0] = float(bp_inp[k]['uv'][0])
              newbp['uv'][1,1] = float(bp_inp[k]['uv'][1])
              newbp['uv'][2,2] = float(bp_inp[k]['uv'][2])
            
        except:        
          pass
       
        try: 
          newbp['expansion'] = float(bp_inp[k]['expansion'])
        except:
          pass         
        try: 
          newbp['rcut'] = units.convert(bp_length, 'ang', float(bp_inp[k]['rcut']))
        except:
          pass      
        
        try:
          newbp['amu_per_crystal'] = float(bp_inp[k]['amu_per_crystal'])
        except:
          pass  
        
        try:
          newbp['b0'] = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['b0']))
        except:
          pass  
        try:
          newbp['e0'] = units.convert(bp_energy, 'EV', bp_inp[k]['e0'])
        except:
          pass        
        try:
          if('ec' in bp_inp[k].keys()):
            newbp['ec'] = numpy.zeros((6,6,),)
            
# CUBIC  C11 C12 C44
            if(len(bp_inp[k]['ec']) == 3):
              c11 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][0]))
              c12 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][1]))
              c44 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][2]))
                            
              newbp['ec'][0,0] = c11
              newbp['ec'][1,1] = c11
              newbp['ec'][2,2] = c11
              
              newbp['ec'][0,1] = c12
              newbp['ec'][0,2] = c12
              newbp['ec'][1,2] = c12
              newbp['ec'][1,0] = c12
              newbp['ec'][2,0] = c12
              newbp['ec'][2,1] = c12
              
              newbp['ec'][3,3] = c44
              newbp['ec'][4,4] = c44
              newbp['ec'][5,5] = c44
              
# ORTHORHOMBIC  C11 C22 C33 C44 C55 C66 C12 C13 C23
            if(len(bp_inp[k]['ec']) == 9):
              c11 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][0]))
              c22 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][1]))
              c33 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][2]))
              c44 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][3]))
              c55 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][4]))
              c66 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][5]))
              c12 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][6]))
              c13 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][7]))
              c23 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][8]))              
              
              newbp['ec'][0,0] = c11
              newbp['ec'][1,1] = c22
              newbp['ec'][2,2] = c33
              
              newbp['ec'][3,3] = c44
              newbp['ec'][4,4] = c55
              newbp['ec'][5,5] = c66
              
              newbp['ec'][0,1] = c12
              newbp['ec'][0,2] = c13
              newbp['ec'][1,2] = c23
              newbp['ec'][1,0] = c12
              newbp['ec'][2,0] = c13
              newbp['ec'][2,1] = c23
              
        except:
          pass
#print(newbp)
# Save to list
        g.bulk_properties.append(newbp)

  def make(label_id, label_str):
# b0 bulk modulus
# e0 cohesive energy  (maybe change to ecoh)
# ec elastic constants
#
   
    bp_d ={
           'bp_key_str': '',
           'label_id': label_id,
           'label_str': label_str,
           'alat': None,
           'uv': numpy.zeros((3,3,),),
           'type': 1,
           'expansion': 4,
           'rcut': 6.5,
           'b0': None,
           'e0': None,
           'ec': None,
           'g': None,
           'e': None,
           'poisson': None,
           'amu_per_crystal': None,
           'type': -1,
           'type_text': '',
          }
    bp_d['uv'][:,:] = 0.0
    bp_d['uv'][0,0] = 1.0
    bp_d['uv'][1,1] = 1.0
    bp_d['uv'][2,2] = 1.0
    return bp_d
           
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def bp_add():  
  
#print(g.pot_labels)
  
    for bp_i in range(len(g.bulk_properties)):
      bp_n = g.bulk_properties[bp_i]
      label_str = bp_n['label_str'].upper()
      if(label_str in g.pot_labels):

        rcut = bp_n['rcut']
        alat = bp_n['alat']
        uv = bp_n['uv']
        label = bp_n['label_id']
        type = bp_n['type']
        expansion = bp_n['expansion']
      
# Add Config
        bp_id = int(bp.add_bp_config(rcut, alat, uv, label, type, expansion))
        g.bulk_properties[bp_i]['fortran_id'] = bp_id

# Add known data
        bp.add_alat(bp_id, bp_n['alat'])
        bp.add_e0(bp_id, bp_n['e0'])
        bp.add_b0(bp_id, bp_n['b0'])
        bp.add_ec(bp_id, bp_n['ec'])
        bp.add_amu_per_crystal(bp_id, bp_n['amu_per_crystal'])
      
        g.bp_ids[bp_id] = {}
        g.bp_ids[bp_id]['label_str'] = label_str
        g.bp_ids[bp_id]['rcut'] = rcut
        g.bp_ids[bp_id]['alat'] = alat
        g.bp_ids[bp_id]['uv'] = uv
        g.bp_ids[bp_id]['label'] = label
        g.bp_ids[bp_id]['type'] = type
        g.bp_ids[bp_id]['expansion'] = expansion
        g.bp_ids[bp_id]['e0'] = bp_n['e0']
        g.bp_ids[bp_id]['b0'] = bp_n['b0']
        g.bp_ids[bp_id]['ec'] = bp_n['ec']
        g.bp_ids[bp_id]['amu_per_crystal'] = bp_n['amu_per_crystal']
        g.bp_ids[bp_id]['type'] = bp_n['type']
        g.bp_ids[bp_id]['type_text'] = bp_n['type_text']
        g.bp_ids[bp_id]['bp_key_str'] = bp_n['bp_key_str']
 
# Add rss multiplication values
    try:
      bp.set_rss_alat(g.inp['rss']['alat'])
    except:
      pass
    try:
      bp.set_rss_e0(g.inp['rss']['e0'])
    except:
      pass
    try:
      bp.set_rss_b0(g.inp['rss']['b0'])
    except:
      pass
    try:
      bp.set_rss_ec(g.inp['rss']['ec'])
    except:
      pass
    try:
      bp.set_rss_g(g.inp['rss']['g'])
    except:
      pass
    try:
      bp.set_rss_e(g.inp['rss']['e'])
    except:
      pass
    try:
      bp.set_rss_v(g.inp['rss']['v'])
    except:
      pass

  @staticmethod
  def bp_output():  
  
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')

      t_pad = 30
      f_pad = 18

      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'LABEL  ' + str(g.bp_ids[bp_id]['label_str'] ), 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Known Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'All units are in ev/Ang unless specified', 1, '', 1)
      std.write_file_line(fh, 'Energy: eV', 1, '', 1)
      std.write_file_line(fh, 'Length: ang', 1, '', 1)
      std.write_file_line(fh, 'Force: eV/ang', 1, '', 1)
      std.write_file_line(fh, 'Pressure: eV/ang3', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.write_file_line(fh, 'Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Calculated Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.write_file_line(fh, 'v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,i,:], f_pad)   
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Compliance:', t_pad, bp.calc_sc[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance (1/GPA):', t_pad, bp.calc_sc_gpa[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc_gpa[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      c11 = 160.230732254e0 * bp.calc_ec[bp_id,0,0]
      c22 = 160.230732254e0 * bp.calc_ec[bp_id,1,1]
      c33 = 160.230732254e0 * bp.calc_ec[bp_id,2,2]
      c44 = 160.230732254e0 * bp.calc_ec[bp_id,3,3]
      c55 = 160.230732254e0 * bp.calc_ec[bp_id,4,4]
      c66 = 160.230732254e0 * bp.calc_ec[bp_id,5,5]
      c12 = 160.230732254e0 * bp.calc_ec[bp_id,0,1]
      c13 = 160.230732254e0 * bp.calc_ec[bp_id,0,2]
      c23 = 160.230732254e0 * bp.calc_ec[bp_id,1,2]
      
      std.write_file_line(fh, 'Stability', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'C11:', t_pad, c11, f_pad)
      std.write_file_line(fh, 'C11C22 - C12C12:', t_pad, (c11*c22)-(c12*c12), f_pad)
      std.write_file_line(fh, 'C11*C22*C33+2*C12*C13*C23-C11*C23*C23-C33*C12*C12:', t_pad, c11*c22*c33+2*c12*c13*c23-c11*c23*c23-c33*c12*c12, f_pad)
      std.write_file_line(fh, 'C44:', t_pad, c44, f_pad)
      std.write_file_line(fh, 'C55:', t_pad, c55, f_pad)
      std.write_file_line(fh, 'C66:', t_pad, c66, f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
          
      std.write_file_line(fh, 'Bulk Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'b0 reuss:', t_pad, bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight:', t_pad, bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg:', t_pad, bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, 'b0 reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Shear Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'G reuss:', t_pad, bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight:', t_pad, bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg:', t_pad, bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, 'G reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight (GPA):', t_pad, 160.230732254e0 * bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg (GPA):', t_pad, 160.230732254e0 * bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Young Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'E:', t_pad, bp.calc_e[bp_id], f_pad)
      std.write_file_line(fh, 'E vec:', t_pad, bp.calc_e_vec[bp_id, :], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Poisson Ratio', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'v:', t_pad, bp.calc_v[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Temperatures', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Melting (K):', t_pad, bp.calc_melting[bp_id], f_pad)
#std.write_file_line(fh, 'Debye:', t_pad, bp.calc_debye[bp_id], f_pad)   # Check, calc might be wrong
      std.write_file_line(fh, '', 1, '', 1) 
             
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss], f_pad)
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS Weighted', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss_w[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss_w[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss_w[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss_w[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss_w[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss_w[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss_w[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total_w[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss_w], f_pad)
      
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      fh.close()

  @staticmethod
  def bp_output_terminal():  
  
    t_pad = 30
    f_pad = 18
    
    for bp_id in range(bp.bp_configs_count):
      
      bp_rcut = g.bp_ids[bp_id+1]['rcut']
      bp_alat = g.bp_ids[bp_id+1]['alat']
      bp_uv = g.bp_ids[bp_id+1]['uv']
      bp_label = g.bp_ids[bp_id+1]['label']
      bp_type = g.bp_ids[bp_id+1]['type']
      bp_expansion = g.bp_ids[bp_id+1]['expansion']
      bp_e0 = g.bp_ids[bp_id+1]['e0']
      bp_b0 = g.bp_ids[bp_id+1]['b0']
      bp_ec = g.bp_ids[bp_id+1]['ec']
      bp_amu_per_crystal = g.bp_ids[bp_id+1]['amu_per_crystal']

      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)   
      std.print_file_line('############################################################################################################', 1, '', 1)
      std.print_file_line('BP KEY          ' + str(g.bp_ids[bp_id+1]['bp_key_str'].upper()), 1, '', 1)
      std.print_file_line('LABEL           ' + str(g.bp_ids[bp_id+1]['label_str']), 1, '', 1)
      std.print_file_line('TYPE            ' + str(g.bp_ids[bp_id+1]['type_text']), 1, '', 1)
      std.print_file_line('############################################################################################################', 1, '', 1)

      std.print_file_line('', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Known Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('All units are in ev/Ang unless specified', 1, '', 1)
      std.print_file_line('Energy: eV', 1, '', 1)
      std.print_file_line('Length: ang', 1, '', 1)
      std.print_file_line('Force: eV/ang', 1, '', 1)
      std.print_file_line('Pressure: eV/ang3', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
      std.print_file_line('Label (element):', t_pad, bp_label, f_pad)
      std.print_file_line('Structure type:', t_pad, bp_type, f_pad)
      std.print_file_line('Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.print_file_line('Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.print_file_line('AMU per crystal:', t_pad, bp_amu_per_crystal, f_pad)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Calculated Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.print_file_line('v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)  
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
    """ 
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')

      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Known Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('All units are in ev/Ang unless specified', 1, '', 1)
      std.print_file_line('Energy: eV', 1, '', 1)
      std.print_file_line('Length: ang', 1, '', 1)
      std.print_file_line('Force: eV/ang', 1, '', 1)
      std.print_file_line('Pressure: eV/ang3', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.print_file_line('Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
      print()
      print()

      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Calculated Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.print_file_line('v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)  
      
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
  """      
        
  @staticmethod
  def bp_eos_plot(dir):       
  
    for bp_id in range(bp.bp_configs_count):
    
# EQUATION OF STATE
    
      s = bp.calc_sizes[bp_id, 0]
    
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(1, 1, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Equation of State')  
      
      plt.xlabel('Volume (ang3)')    
      plt.ylabel('Energy (eV)')
          
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies[bp_id, 0, 0:s], color='k',  marker="x", ls='')
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies_fit[bp_id, 0, 0:s], color='k', ls='solid')

      plt.savefig(dir + '/' + 'eos_' + str(bp_id) + '.svg')
      plt.savefig(dir + '/' + 'eos_' + str(bp_id) + '.eps')
      
# ELASTIC CONSTANTS
      
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(3, 3, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Elastic Constant Curves')    
    
      for dn in range(9):
        s = bp.calc_sizes[bp_id, dn + 1]
    
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies[bp_id, dn + 1, 0:s],
                                                 color='k',  marker="x", ls='')
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies_fit[bp_id, dn + 1, 0:s],
                                                 color='k', ls='solid')
        axs[int(numpy.floor(dn/3)), dn % 3].set_title('Distortion D' + str(dn + 1))
        axs[int(numpy.floor(dn/3)), dn % 3].set_xlabel('Strain (Expanded Alat)')
        axs[int(numpy.floor(dn/3)), dn % 3].set_ylabel('Energy (eV)')
               
      plt.savefig(dir + '/' + 'ec_' + str(bp_id) + '.svg')
      plt.savefig(dir + '/' + 'ec_' + str(bp_id) + '.eps')
      
    """
    n = globals.d['eos_data_size']
    x = np.linspace(globals.d['eos_data'][1,0], globals.d['eos_data'][1,n-1], 101)    
    y = np.zeros((101,),)
    p = numpy.zeros((4,),)
    p[:] = globals.d['eos_fitting'][:]
    y[:] = run_eos.bm_calc(x[:], p)
    
    n = globals.d['eos_data_size'] 
    plt.plot(globals.d['eos_data'][1,:n], globals.d['eos_data'][2,:n], color='k',  marker="x", ls='')
    plt.plot(x[:], y[:], color='k', ls='solid')
    """
  
###########################################
#  CLASS unit
###########################################
class units:
  
  @staticmethod
  def convert(conv_from, conv_to, value_in):
    try:
      value_in = float(value_in)
    except:
      return None
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# LENGTH METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }
    
# AREA METERS SQUARED
    area = {
    'M2': 1.0,
    'CM2': 1E4,
    'MM2': 1E6,
    'UM2': 1E12,
    'NM2': 1E18,
    'ANG2': 1E20,
    }
    
# VOLUME METERS CUBED
    volume = {
    'M2': 1.0,
    'CM2': 1E6,
    'MM2': 1E9,
    'UM2': 1E18,
    'NM2': 1E27,
    'ANG2': 1E30,
    }

# ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    }

# FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
# VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
# PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
# CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }
    
# TEMPERATURE
    
    unit_list = [length, area, volume, energy, force, velocity, pressure, charge_density]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  
"""  
  @staticmethod
  def convert(conv_from, conv_to, value_in):

    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

# METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

# J
    energy = {
    'J': 1.0,
    'EV': 6.242E18,
    'RY': 4.5874E17,
    }

    if(conv_from in length.keys() and conv_to in length.keys()):
      return round((length[conv_to] / length[conv_from]) * float(value_in),9)

    if(conv_from in energy.keys() and conv_to in energy.keys()):
      return round((energy[conv_to] / energy[conv_from]) * float(value_in),9)
"""

###########################################
#  CLASS efs_cal
###########################################
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

###########################################
#  CLASS bp_cal
###########################################
class bp_calc:

  def run():  
  
    print("Calc Bulk Properties") 

    bp_calc.output_dir()
    
# Setup BP
    bp_calc.init()
    
# Load potentials
    potential.bp_add_potentials()
    
# Add required BP structures (add, make ghost configs, make NLs)
    b_props.bp_add()
        
# Calculate energies of configurations for all structures
    bp.energy() 
   
#print(bp.cc)
#print(bp.calc_count)
#print(bp.cc_log[0,1)

# Calculate BP
    bp.calculate_bp()   
    
#print(bp.rss)
    
# Output to Terminal
    b_props.bp_output_terminal()
    
# Output to File
#b_props.bp_output()

    b_props.bp_eos_plot(g.dirs['wd'] + '/bp')
    
# Plots
#b_props.bp_eos_plot()
#potential.plot_fortran_potentials()
#potential.plot_python_potentials()
    
  def output_dir():
    std.make_dir(g.dirs['wd'] + '/bp')
    
  def init():
# Log
    main.log('BP Allocated Memory\n')
    try:
      mem = str(g.inp['mem']['bp'])
    except:
      mem = "500MB"
    main.log(mem + '\n')
    mem = std.mem_value(mem)    
    g.memory['bp']['c'] = int(1 * (mem / 7840))
    g.memory['bp']['g'] = int(12 * (mem / 7840))
    g.memory['bp']['nl'] = int(100 * (mem / 7840))
    main.log(str(g.memory['bp']['c']) + " " + str(g.memory['bp']['g']) + " " + str(g.memory['bp']['nl']) + '\n')
    
# Initialise
    bp.init(g.memory['bp']['c'], g.memory['bp']['g'], g.memory['bp']['nl'])
 
    bp.eos_size = g.bp_input['eos_size']
    bp.eos_strain = g.bp_input['eos_strain']
    bp.ec_size = g.bp_input['ec_size']
    bp.ec_strain = g.bp_input['ec_strain']
    
    """g.bp_input = {
    'dir': '',
    'bp_file': None,
    'eos_size': 10,
    'eos_strain': 10,
    'ec_size': 10,
    'ec_strain': 10,
    }"""
    
  def set_weights():
    bp.set_rss_total(g.rss_weights['bp'])
    bp.set_rss_a0(g.rss_weights['a0'])
    bp.set_rss_e0(g.rss_weights['e0'])
    bp.set_rss_b0(g.rss_weights['b0'])
    bp.set_rss_ec(g.rss_weights['ec'])
    bp.set_rss_g(g.rss_weights['g'])
    bp.set_rss_e(g.rss_weights['e'])
    bp.set_rss_v(g.rss_weights['v'])
    bp.set_rss_neg_ec(g.rss_weights['negec'])
    
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

    """
# Bias so worse RSS if the density is out of range
    f = 1.0  
    if(bp.max_density < g.rss_max_density['min']):
      f = 1.0 + (100.0 * (g.rss_max_density['min'] - bp.max_density))**4
    if(bp.max_density > g.rss_max_density['max']):
      f = 1.0 + (g.rss_max_density['scale_factor'] * (bp.max_density - 0.8))**g.rss_max_density['scale_exponent']
    """
    f = 1.0  
    if(bp.max_density == 0.0):
      f = g.rss_max_density['zero_density_factor']  

    try:
# Make Dictionary
      g.rss['bp'] = {}
      g.rss['bp']['ok'] = True  
      g.rss['bp']['cc'] = int(bp.cc)
# Totals
      g.rss['bp']['total_rss'] = float(bp.rss_total_rss)
      g.rss['bp']['total_rss_weighted'] = f * float(bp.rss_total_rss_w)
# Individual RSS
      g.rss['bp']['a0'] = float(bp.rss_by_type[0])
      g.rss['bp']['e0'] = float(bp.rss_by_type[1])
      g.rss['bp']['b0'] = float(bp.rss_by_type[2])
      g.rss['bp']['ec'] = float(bp.rss_by_type[3])
      g.rss['bp']['g'] = float(bp.rss_by_type[4])
      g.rss['bp']['e'] = float(bp.rss_by_type[5])    
      g.rss['bp']['v'] = float(bp.rss_by_type[6])
# Weighted RSS
      g.rss['bp']['a0_weighted'] = f * float(bp.rss_by_type_w[0])
      g.rss['bp']['e0_weighted'] = f * float(bp.rss_by_type_w[1])
      g.rss['bp']['b0_weighted'] = f * float(bp.rss_by_type_w[2])
      g.rss['bp']['ec_weighted'] = f * float(bp.rss_by_type_w[3])
      g.rss['bp']['g_weighted'] = f * float(bp.rss_by_type_w[4])
      g.rss['bp']['e_weighted'] = f * float(bp.rss_by_type_w[5])    
      g.rss['bp']['v_weighted'] = f * float(bp.rss_by_type_w[6])
    except:
# Make Dictionary
      g.rss['bp'] = {}
      g.rss['bp']['ok'] = False
      g.rss['bp']['cc'] = 0 
           
    for i in range(bp.residuals_n):
      g.rss['residual'].append(f * float(bp.residuals[i]))
      
###########################################
#  CLASS es_cal
###########################################
class es_calc:

  def run():  
  
    print("Energies: Surface, Vacany etc") 

# Setup ES
    es.init()
    potential.es_add_potentials()
    es_calc.es_add()
    
#es.energy()
#for i in range(es.cc):
#  print(es.config_energy[i,:])
    es.calculate_es()
    
  @staticmethod
  def es_add():  
  
    print(g.es)
    for item in g.es:
      es_id = es.add_es_config(item['rcut'], item['alat'], item['label_id'], item['type'])
  
    """ 
# Test values
    
    rcut = 6.5
    alat = 4.04
    label = 1
    type = 3
      
# Add Config
    es_id = es.add_es_config(rcut, alat, label, type)
#print(es_id)
    
    rcut = 6.5
    alat = 3.5
    label = 1
    type = 2
    es_id = es.add_es_config(rcut, alat, label, type)
#print(es_id)
    """

  @staticmethod
  def load():
    if('es' not in g.inp.keys()):
      return None
    if('es_file' not in g.inp['es'].keys()):
      return None
    try: 
      dir = g.inp['es']['dir'].strip()
      es_file = std.path(g.inp['es']['dir'], g.inp['es']['es_file'])
    except:
      dir = ""
      es_file = g.inp['es']['es_file']
    
# Read BP data file
    es_inp = read_config.read_file(es_file)

# Make list
    g.es = []
    
# READ IN UNITS
    try:
      es_pressure = es_inp['units']['pressure']
    except:
      es_pressure = 'GPA'
    try:
      es_length = es_inp['units']['length']
    except:
      es_length = 'ang'
    try:
      es_energy = es_inp['units']['energy']
    except:
      es_energy = 'ev'
      
    for k in es_inp.keys():
      if('potlabel' in es_inp[k].keys() and 'alat' in es_inp[k].keys()):
        potlabel = es_inp[k]['potlabel']
        label_str, label_id = labels.add(potlabel)     
        
        new_es = es_calc.make(label_id, label_str)

# ALAT
# There must be an alat value set
        new_es['alat'] = units.convert(es_length, 'ang', float(es_inp[k]['alat']))        
          
# UNIT VECTOR
        new_es['uv'][:,:] = 0.0
        new_es['uv'][0,0] = 1.0
        new_es['uv'][1,1] = 1.0
        new_es['uv'][2,2] = 1.0        
        try: 
          if('uv' in es_inp[k].keys()):
# CUBIC
            if(len(es_inp[k]['uv']) == 1):
              new_es['uv'][0,0] = float(es_inp[k]['uv'][0])
              new_es['uv'][1,1] = float(es_inp[k]['uv'][0])
              new_es['uv'][2,2] = float(es_inp[k]['uv'][0])
# CUBIC
            if(len(es_inp[k]['uv']) == 3):
              new_es['uv'][0,0] = float(es_inp[k]['uv'][0])
              new_es['uv'][1,1] = float(es_inp[k]['uv'][1])
              new_es['uv'][2,2] = float(es_inp[k]['uv'][2])            
        except:        
          pass
          
# TYPE
        try:
          if(es_inp[k]['type'].lower() == 'sc'):
            new_es['type'] = 1
          elif(es_inp[k]['type'].lower() == 'bcc'):
            new_es['type'] = 2
          elif(es_inp[k]['type'].lower() == 'fcc'):
            new_es['type'] = 3
          elif(es_inp[k]['type'].lower() == 'zb'):
            new_es['type'] = 4
        except:
          pass 
          
# EXPANSION
        try: 
          new_es['expansion'] = float(es_inp[k]['expansion'])
        except:
          pass   
        
# RCUT
        try: 
          new_es['rcut'] = units.convert(es_length, 'ang', float(es_inp[k]['rcut']))
        except:
          pass  

# RCUT
        try: 
          new_es['surface_energy'] = units.convert(es_length, 'ang', float(es_inp[k]['surface_energy']))
        except:
          pass 
       
# Add to list
        g.es.append(new_es)  

#
      
  def make(label_id, label_str):
# b0 bulk modulus
# e0 cohesive energy  (maybe change to ecoh)
# ec elastic constants
#
   
    es_d ={
           'label_id': label_id,
           'label_str': label_str,
           'alat': None,
           'uv': numpy.zeros((3,3,),),
           'type': 1,
           'expansion': 4,
           'rcut': 6.5,
           'surface_energy': None,
          }
    es_d['uv'][:,:] = 0.0
    es_d['uv'][0,0] = 1.0
    es_d['uv'][1,1] = 1.0
    es_d['uv'][2,2] = 1.0
    return es_d
               
###########################################
#  CLASS rss_cal
###########################################
class rss_calc:

  def run():
    print("Calc RSS") 

# Make dir
    rss_calc.output_dir()   

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
 
      print('BP  ' + str(bp_id)) 
      print('================') 

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
    print('RSS: ' + str(rss))
    print('')
    print('')
    print('Max Density:', bp.max_density, efs.max_density)
    print('')
    print('')
#print(g.rss)
    print('')
    
    exit()

# Plots
#potential.plot_fortran_potentials()
#potential.plot_python_potentials()
    potential.plot_python_potentials(g.dirs['wd'] + '/rss/plots/potential_python')
    potential.plot_fortran_potentials(g.dirs['wd'] + '/rss/plots/potential_fortran')
    
    b_props.bp_output()
    b_props.bp_eos_plot(g.dirs['wd'] + '/rss/plots')

  def output_dir():
    std.make_dir(g.dirs['wd'] + '/rss')  
    std.make_dir(g.dirs['wd'] + '/rss/plots')
    
# ASSUMES EFS AND BP ALREADY SET UP
  def run_calc(save_in_top=True): 
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
    g.rss['current_bp_max_density'] = copy.deepcopy(bp.max_density)
    g.rss['current_efs_max_density'] = copy.deepcopy(efs.max_density)
    if(g.rss['best'] == None or g.rss['current'] < g.rss['best']):
      main.log('Best rss: ' + str(g.rss['best']) + ' (since last: ' + str(g.rss['since_improvement']) + ')')
      g.rss['best_bp_max_density'] = copy.deepcopy(bp.max_density)
      g.rss['best_efs_max_density'] = copy.deepcopy(efs.max_density)
      g.rss['best'] = g.rss['current']
      g.rss['since_improvement'] = 0
      
      g.efs_results_best = copy.deepcopy(g.efs_results)
      g.bp_results_best = copy.deepcopy(g.bp_results)
    
    """
# TOP 100 List
    if(save_in_top):
      list_len = g.fitting['top_parameters']

# COPY ARRAYS
      p_now = numpy.copy(potential.get_parameters())
      efs_results = copy.deepcopy(g.efs_results)
      bp_results = copy.deepcopy(g.bp_results)

      if(len(g.top_parameters) == 0):
        g.top_parameters.append([rss, p_now, efs_results, bp_results, bp.max_density, efs.max_density])
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
          g.top_parameters.insert(i, [rss, p_now, efs_results,bp_results, bp.max_density, efs.max_density])
        if(len(g.top_parameters) > list_len):
          g.top_parameters.pop()
    """ 
      
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

###########################################
#  CLASS p
###########################################
class pf:

  start_time = 0.0
  end_time = 0.0  
  
  def run():   
  
# Start Time
    pf.start_time = time.time()
  
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
    p = potential.get_start_parameters()
    pf_potential.update(p)
    rss = pf.get_rss()
    pf.set_current(p)

# Load saved
    pf_saved.run()

#pf_saa.run()

    for fit in g.fit:
      pf.fit = fit
      if(pf.fit['type'].lower() == "random"):        
        pf_random.run()
      elif(pf.fit['type'].lower() == "sa"): 
        pf_sa.run()
      elif(pf.fit['type'].lower() == "ga"): 
        pf_genetic.run()
      elif(pf.fit['type'].lower() == "ng"): 
        pf_ng.run()

    pf_top.save(g.dirs['wd'] + '/save', 'top_parameters.txt')

# Output
    potential_output.full()
          
# End Time
    pf.end_time = time.time()
  
# Log fit time
    main.log("Fit time: ", str(pf.end_time - pf.start_time))

    print(g.top_parameters[0][0])
    print(g.top_parameters[0][1])
    print(g.top_parameters[0][2])
    print(g.top_parameters[0][3])
    print(g.top_parameters[0][4])
    print(g.top_parameters[0][5])

# Final
    pf_final.run()

#pf_ng.run()

# Enhance top
#pf_setop.run()

# Run Genetic
#pf_genetic.run()
    
# Run Sim Annealing on Best Result
#pf_sa.run()

#pf_gd.run()
#pf_steps.run()
  
#pf_ng.run()

# Display
#g.pfdata['stage'] = 'Finished'
#display.finish()
      
######################################################
# SET UP - initialise EFS, BP and any other modules
#          and run first calc
######################################################

  def set_current(p):
    g.pfdata['p']['current'] = numpy.copy(p)
    
  def get_rss(save_in_top=True):  
    rss = rss_calc.run_calc(save_in_top)  
        
    g.pfdata['rss']['current'] = rss
    g.pfdata['rss']['counter'] += 1
    g.pfdata['max_density']['bp_current'] = bp.max_density
    g.pfdata['max_density']['efs_current'] = efs.max_density
    
# If successful
    if(rss is not None):   
    
# Store Details
      g.pfdata['rss']['counter_successful'] += 1
      g.pfdata['rss']['since_improvement'] += 1    
      if(g.pfdata['rss']['start'] == None):
        g.pfdata['rss']['start'] = rss
      if(g.pfdata['rss']['best'] == None or rss < g.pfdata['rss']['best']):
        g.pfdata['p']['best'] = numpy.copy(g.pfdata['p']['current'])
        g.pfdata['rss']['best'] = rss
        g.pfdata['bp']['best'] = copy.deepcopy(g.bp_results)        
        g.pfdata['rss']['since_improvement'] = 0
        g.pfdata['max_density']['bp_best'] = bp.max_density
        g.pfdata['max_density']['efs_best'] = efs.max_density
        
# Save in top
      if(save_in_top):
        list_len = g.fitting['top_parameters']

# COPY ARRAYS
        p_now = numpy.copy(potential.get_parameters())
        efs_results = copy.deepcopy(g.efs_results)
        bp_results = copy.deepcopy(g.bp_results)
        rss_details = copy.deepcopy(g.rss)

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
    pf_display.output()    
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
    
###########################################
#  CLASS g
###########################################
class gd:

  precision = 1.0e-6
  max_iterations = 10
  h = 1.0e-8
  momentum = 0.1
  p_in = None
  rss_in = None
  p_out = None
  rss_out = None
  fix = None

# Central Difference
  def df(p):
    d = numpy.zeros((len(p),),)
    for i in range(len(p)):
      if(gd.fix[i] == 0.0):
        d[i] = 0.0
      else:
        p_f = numpy.copy(p)
        p_f[i] = p_f[i] + gd.h
        p_b = numpy.copy(p)
        p_b[i] = p_b[i] - gd.h
        d[i] = (gd.rss(p_f) - gd.rss(p_b)) / (2 * gd.h)
    return d

  def line_search(p, dp):
# Back track
    best_rss = gd.rss_in
    best_gamma = 0.0
    gamma = 20.0 * gd.last_gamma 
    if(gamma <1.0e-8):
      gamma <1.0e-6
    loop = True
    n = 0
    while(gamma > 1.0e-12):
      n = n + 1
      p_test = p - gamma * dp
      rss = gd.rss(p_test)
      gamma = 0.5 * gamma
      if(best_rss is None or rss < best_rss):
        p_best = p_test
        best_rss = rss
        best_gamma = gamma
    gd.last_gamma = best_gamma
    return best_rss, best_gamma * dp
    
# Gradient Descent
  def opt(f_rss, p0, fix=None):
    gd.rss = f_rss   
    gd.last_gamma = 1.0
    p = p0
    
    if(type(fix) == numpy.ndarray):
      if(len(p) == len(fix)):
        for i in range(len(fix)):
          if(fix[i] != 0.0):
            fix[i] = 1.0
      else:
        fix = numpy.zeros((len(p0),),)
        fix[:] = 1.0
    else:
      fix = numpy.zeros((len(p0),),)
      fix[:] = 1.0
    gd.fix = fix  
    gd.p_in = numpy.copy(p)
    gd.rss_in = gd.rss(p)
    
    best_p = numpy.copy(gd.p_in)
    best_rss = gd.rss_in

    n = 0
    last_dp = 0.0
    while(n < gd.max_iterations):
      df = gd.df(p)
      rss, dp = gd.line_search(p, df)
      if(rss > best_rss):
        n = gd.max_iterations
      else:
        p = p - (gd.momentum * last_dp + dp)
        last_dp = dp
        rss = gd.rss(p)
        best_p = numpy.copy(p)
        best_rss = rss
        n = n + 1
        
    gd.p_out = numpy.copy(best_p)
    gd.rss_out = best_rss
    
    return p, rss

###########################################
#  CLASS pf_star
###########################################
class pf_start:

  def run():
  
    pf_display.stage = 'START - LOAD SAVED'
    
    p = potential.get_start_parameters()
    pf_potential.update(p)
    rss = pf.get_rss()

#pf_pgradient.run()

# Load saved parameters
    params = pf_top.load(g.dirs['wd'] + '/save', 'top_parameters.txt')    
    for p in params:
      pf_potential.update(p)
      rss = pf.get_rss()

###########################################
#  CLASS pf_step
###########################################
class pf_steps:

  def run():

    p = numpy.copy(g.pfdata['params']['start'])
    rss = pf.get_rss()
    
    for n in range(5):
      p_new = numpy.copy(p)
      p_new = p_new + pf_steps.random_step()

      pf_potential.update(p_new)
      rss_new = pf.get_rss()
  
      if(rss_new < rss):
        p = numpy.copy(p_new)
        rss_new = rss

  def random_step(m = 0.01):
# Get random parameters 0 to 1
    return m * (0.5 - numpy.random.rand(g.pfdata['params']['count']))
    
  def process_list(p_list, steps = 100):

    p_new_list = []

    for n in range(len(p_list)):

      g.pfdata['stage'] = 'STEP_' + str(n+1)

      p_current = numpy.copy(p_list[n])
      rss_current = pf.get_rss()
      rss_start = rss_current

      for n in range(steps):
        p_new = numpy.copy(p_current)
        p_new = p_new + pf_steps.random_step(0.0001)

        pf_potential.update(p_new)
        rss = pf.get_rss(False)

        if(rss < rss_current):
          p_current = p_new
          rss_current = rss
      
      p_new_list.append(p_current)

    return p_new_list

###########################################
#  CLASS pf_geneti
###########################################
class pf_genetic:

  count = 0
  start_time = 0
  variation_factor = 1.0

  def run():    
    if(pf.fit['gens'] == 0):
      return 0
    
    g.pfdata['stage'] = 'Genetic Algorithm ' + str(pf_genetic.count)
    g.pfdata['stage_brief'] = 'GA' + str(pf_genetic.count)

#
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 
    pf_genetic.start_time = time.time()
    pf_genetic.count = pf_genetic.count + 1
  
    main.log_title("Genetic Fit")   
    
# Load from input
    pf_genetic.width = g.pfdata['psize']
    pf_genetic.pop_size = pf.fit['pop_size']
    pf_genetic.fresh_size = pf.fit['fresh_size']
    pf_genetic.generations = pf.fit['gens']
    pf_genetic.no_clone_var = pf.fit['no_clone_var']
    pf_genetic.gen_variation_multiplier = pf.fit['gen_variation_multiplier']

# Force to be even
    if(pf_genetic.pop_size % 2 != 0):
      pf_genetic.pop_size = pf_genetic.pop_size + 1
    if(pf_genetic.fresh_size % 2 != 0):
      pf_genetic.fresh_size = pf_genetic.fresh_size + 1

# Calculate other sizes
    pf_genetic.children_size = pf_genetic.pop_size + 2 * pf_genetic.fresh_size
    pf_genetic.merge_size = 2 * pf_genetic.pop_size + 3 * pf_genetic.fresh_size
    
# Set Up Pop + Fresh arrays
    pf_genetic.pop = numpy.zeros((pf_genetic.pop_size, pf_genetic.width,),)
    pf_genetic.pop_rss = numpy.zeros((pf_genetic.pop_size,),)
    pf_genetic.fresh = numpy.zeros((pf_genetic.fresh_size, pf_genetic.width,),)
    pf_genetic.fresh_rss = numpy.zeros((pf_genetic.fresh_size,),)
    pf_genetic.children = numpy.zeros((pf_genetic.children_size, pf_genetic.width,),)
    pf_genetic.children_rss = numpy.zeros((pf_genetic.children_size,),)
    pf_genetic.pop_merged = numpy.zeros((pf_genetic.merge_size, pf_genetic.width,),)
    pf_genetic.pop_merged_rss = numpy.zeros((pf_genetic.merge_size,),)
  
# Parents array
    pf_genetic.parents = numpy.arange(pf_genetic.pop_size)

# Initialise population
    g.pfdata['stage'] = 'Genetic Fit ' + str(pf_genetic.count) + ' Initialising'
  
# Copy out parameters from "top_parameters"
    start_parameters = []
    if(len(g.top_parameters) > 0):
      for i in range(len(g.top_parameters)):
        new_p = numpy.copy(g.top_parameters[i][1][:])
        start_parameters.append(new_p)
    else:
      p = potential.get_start_parameters()
      new_p = numpy.copy(p)
      start_parameters.append(new_p)

# Load start parameters into pop, and fill remaining with random parameters
    n = 0
    for pn in range(pf_genetic.pop_size):
      loop = True
      while(loop): 
        n = n + 1
        if(n <= len(start_parameters)):
          pf_genetic.pop[pn, :] = start_parameters[n-1][:]
        else:
          pf_genetic.pop[pn, :] = pf_parameters.random_p()
        
# Try - if it fails or rss == None, try next
        pf_potential.update(pf_genetic.pop[pn, :])
        try:
# Update
          rss = pf.get_rss()
          if(rss is not None):
            loop = False
            pf_genetic.pop_rss[pn] = rss
        except:
          pass

###################################
# LOOP THROUGH GENERATIONS
###################################

    pf_genetic.generation = 0
    for gen in range(pf_genetic.generations):
      pf_genetic.generation = pf_genetic.generation + 1
      pf_genetic.run_gen()

# End
    pf.summary_line("GENETIC ALGORITHM " + str(pf_genetic.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("GENETIC_ALGORITHM_" + str(pf_sa.count))
  
  def run_gen():

    ss = 'Genetic Fit ' + str(pf_genetic.count) + ' Gen ' + str(pf_genetic.generation) + ' '  

# Shuffle parents array
    numpy.random.shuffle(pf_genetic.parents) 

    c = 0
# Breed Population
    g.pfdata['stage'] = ss + 'Breed Population'
    for p in range(pf_genetic.pop_size // 2):  
      loop = True
      while(loop):
        loop = pf_genetic.breed_event(p, c, 'p+p')      
      c = c + 2

# Make Fresh
    g.pfdata['stage'] = ss + 'Initialising Fresh Population'   
    pf_genetic.make_fresh()

# Shuffle parents array
    numpy.random.shuffle(pf_genetic.parents) 

# Breed Population-Fresh
    g.pfdata['stage'] = ss + 'Breed With Fresh Population'  
    for p in range(pf_genetic.fresh_size):  
      loop = True
      while(loop):
        loop = pf_genetic.breed_event(p, c, 'p+f')      
      c = c + 2

# Merge populations and select top to form next generation
    g.pfdata['stage'] = ss + 'Merge Populations'   
    pf_genetic.merge()

# Load best back into population
    pf_genetic.pop[:, :] = pf_genetic.pop_merged[0:pf_genetic.pop_size, :]
    pf_genetic.pop_rss[:] = pf_genetic.pop_merged_rss[0:pf_genetic.pop_size]
   
# Change variation for next generation
    pf_genetic.variation_factor = pf_genetic.variation_factor * pf_genetic.gen_variation_multiplier

  def breed_event(p, c, opt):
    pc = pf_genetic.width
    
    pa = pf_genetic.parents[p]
    if(opt == 'p+p'):  
      pb = pf_genetic.parents[p + pf_genetic.pop_size // 2]
      pb_array = 'pop'
    if(opt == 'p+f'): 
      pb = p
      pb_array = 'fresh'
    
    ca = c
    cb = c + 1
    
# Breed
    state = pf_genetic.get_state()
    for i in range(pc):
      state = pf_genetic.get_state(state)
      if(state):
        pf_genetic.children[ca, i] = pf_genetic.pop[pa, i]
        if(pb_array == 'pop'):
          pf_genetic.children[cb, i] = pf_genetic.pop[pb, i]
        elif(pb_array == 'fresh'):
          pf_genetic.children[cb, i] = pf_genetic.fresh[pb, i]
      else:     
        if(pb_array == 'pop'):
          pf_genetic.children[ca, i] = pf_genetic.pop[pb, i]
        elif(pb_array == 'fresh'):
          pf_genetic.children[ca, i] = pf_genetic.fresh[pb, i]
        pf_genetic.children[cb, i] = pf_genetic.pop[pa, i]
      
# Mutate
#g.pfdata['params']['children'][ca, :] = pf_genetic.mutate(g.pfdata['params']['children'][ca, :], g.fit['mutate_chance'])
#g.pfdata['params']['children'][cb, :] = pf_genetic.mutate(g.pfdata['params']['children'][cb, :], g.fit['mutate_chance'])
    
# No clones - Child A
    for pn in range(pf_genetic.pop_size):
      if((pf_genetic.children[ca, :].all() == pf_genetic.pop[pn,:]).all()):
        r = numpy.random.rand(pc)   
        pf_genetic.children[ca, :] = pf_parameters(pf_genetic.children[ca, :],  pf_genetic.no_clone_var)
        break
        
# No clones - Child B
    for pn in range(pf_genetic.pop_size):
      if((pf_genetic.children[cb, :].all() == pf_genetic.pop[pn,:]).all()):
        r = numpy.random.rand(pc)   
        pf_genetic.children[cb, :] = pf_parameters(pf_genetic.children[cb, :],  pf_genetic.no_clone_var)
        break

# Check the children actually give valid parameters
    loop = False     
    pf_potential.update(pf_genetic.children[ca, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
    pf_potential.update(pf_genetic.children[cb, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
    if(loop == False):
      pf_genetic.children_rss[ca] = rss
      pf_genetic.children_rss[cb] = rss
      
# Loop again or not
    return loop
    
  def get_state(state = None):
    if(state == None):
      if(random.uniform(0.0, 100.0) > 50.0):
        return True
      return False
    else:
      if(random.uniform(0.0, 100.0) > 50.0):
        if(state):
          return False
        return True
      return state

  def make_fresh():    

    for p in range(pf_genetic.fresh_size):
      loop = True
      while(loop):
        try:
          rn = min(len(g.top_parameters) - 1, numpy.floor((len(g.top_parameters) + 1) * random.uniform(0.0, 1.0)**3))
          pbest = numpy.copy(g.top_parameters[0][1])
          pf_genetic.fresh[p, :] = pf_parameters.random_p(pbest, pf_genetic.variation_factor)

#pf_genetic.fresh[p, :] = numpy.copy(g.top_parameters[0][1])

          pf_potential.update(pf_genetic.fresh[p, :])
          rss = pf.get_rss()
          if(rss is not None):
            loop = False
            pf_genetic.fresh_rss[p] = rss
        except:
          pass 

# Used to merge parents, fresh and all children into ordered array
  def merge():
    ps = pf_genetic.pop_size
    fs = pf_genetic.fresh_size
    cs = pf_genetic.pop_size + 2 * pf_genetic.fresh_size
      
# Reset array
    pf_genetic.pop_merged[:,:] = 0.0
    pf_genetic.pop_merged[0:ps,:] = pf_genetic.pop[:,:]   
    pf_genetic.pop_merged[ps:ps+fs,:] = pf_genetic.fresh[:,:]   
    pf_genetic.pop_merged[ps+fs:ps+fs+cs,:] = pf_genetic.children[:,:]  
    pf_genetic.pop_merged_rss[:] = 0.0
    pf_genetic.pop_merged_rss[0:ps] = pf_genetic.pop_rss[:] 
    pf_genetic.pop_merged_rss[ps:ps+fs] = pf_genetic.fresh_rss[:] 
    pf_genetic.pop_merged_rss[ps+fs:ps+fs+cs] = pf_genetic.children_rss[:] 
        
# Sort
    pf_genetic.sort_merged()
    
  def sort_merged():
    sort.sort_1d_dp_asc(pf_genetic.pop_merged_rss[:])
    pf_genetic.pop_merged_rss = sort.apply_keytable_1d_dp(pf_genetic.pop_merged_rss)
    pf_genetic.pop_merged = sort.apply_keytable_2d_dp(pf_genetic.pop_merged)

###########################################
#  CLASS displa
###########################################
class display:

  display_header_set = False
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
#print(g.rss)
    
  def output_1():
  
# Best BP
    bp_best = g.pfdata['bp_best']
    
    if(display.display_header_set == False or display.last_stage == None or display.last_stage != g.pfdata['stage']):
      display.last_stage = g.pfdata['stage']
      display.display_header_set = True 
      
# PRINT HEADER
      display.print_line()
      print("# Cycles:           " + str(g.pfdata['cycle']['total_cycles']))
      print("# Generations:      " + str(g.pfdata['generation']['total_generations']))
      print("# Pop Size:         " + str(g.fit['pop_size']))
      print("# Fresh Size:       " + str(g.fit['fresh_size']))
      print("# Configs/second:   " + str(g.benchmark['configspersec']))
      print("# Stage:            " + g.pfdata['stage'])
      display.print_line()
  
    print("##### ")
    print("# Cycle/Gen: ", display.pad_r(g.pfdata['cycle']['counter'], 5), display.pad_r(g.pfdata['generation']['counter'], 5), end=" ")
    print("Counter: ", display.pad_r(g.pfdata['rss']['counter'], 16), end=" ")
    print("# Configs/second:   " + display.pad_r(g.benchmark['configspersec'], 16), end=" ")
    print("# Atoms/second:     " + display.pad_r(g.benchmark['atomspersec'], 16), end=" ")
    print()
    print("# Timer:            " + display.pad_l(time.time() - g.times['start'], 16))
    print("# RSS: ", display.pad_r(g.pfdata['rss']['current'], 20), end="          ")
    print("Best: ", display.pad_r(g.pfdata['rss']['best'], 20), end=" ")
    print()    
    print("Max Density: ", display.pad_r(bp.max_density, 23), display.pad_r(efs.max_density, 23), end=" ")
    print()
    print("#        ID   a0       e0       B0       C11      C12      C44 ")
    try:
      for bp_id in g.bp_results['input'].keys():
        print(  display.pad_r("# Calc: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in bp_best['input'].keys():
        print(  display.pad_r("# Best: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in g.bp_known['input'].keys():      
        print(  display.pad_r("# Known: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    
  def output_2():
    if(display.last_stage == None or display.last_stage != g.pfdata['stage']):
      display.last_stage = g.pfdata['stage']
      print("Stage: " + g.pfdata['stage'])
  
    print("   " , g.pfdata['rss']['counter'], g.pfdata['rss']['current'], g.pfdata['rss']['best'])
  
  def output_3(results=None):

# HEADER
    display.header_3()

# PRINT RSS
    print("# RSS:                ", g.pfdata['rss']['current'], "  [",g.pfdata['rss']['best'],"]   (Max density: ", g.rss['current_bp_max_density'],  g.rss['current_bp_max_density'], ")") 
    print("# RSS (BEST):         ", g.pfdata['rss']['best'], "   (Max density: ", g.rss['best_bp_max_density'],  g.rss['best_bp_max_density'], ")") 
    display.print_line()

# Best BP
    bp_best = g.pfdata['bp_best']
    
# PRINT VALUES
    print("#        ID   a0       e0       B0       C11      C12      C44 ")
    try:
      for bp_id in g.bp_results['input'].keys():
        print(  display.pad_r("# Calc: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_results['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in bp_best['input'].keys():
        print(  display.pad_r("# Best: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in g.bp_known['input'].keys():      
        print(  display.pad_r("# Known: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    display.print_line()
    
# PRINT PARAMETERS
    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("Fn:" + str(fn) + "[S] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,0],8), end='')
        print()
        print("         ", end="")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,1],8), end='')
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
    
  def finish_1(results=None):
# HEADER
# PRINT HEADER
    display.print_line()
    print("# Timer:            " + display.pad_l(time.time() - g.times['start'], 16))
    print("# Cycles:           " + str(g.pfdata['cycle']['total_cycles']))
    print("# Generations:      " + str(g.pfdata['generation']['total_generations']))
    print("# Pop Size:         " + str(g.fit['pop_size']))
    print("# Fresh Size:       " + str(g.fit['fresh_size']))
    print("# Configs/second:   " + str(g.benchmark['configspersec']))
    print("# Stage:            " + g.pfdata['stage'])
    display.print_line()
    
# PRINT RSS
    print("# BEST RSS:                  ",g.pfdata['rss']['best'],"") 
    display.print_line()
    
# Best BP
    bp_best = g.pfdata['bp_best']
    
# PRINT VALUES
    print("#        ID   a0       e0       B0       C11      C12      C44 ")
    try:
      for bp_id in bp_best['input'].keys():
        print(  display.pad_r("# Best: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in g.bp_known['input'].keys():      
        print(  display.pad_r("# Known: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    display.print_line()
    
# PRINT PARAMETERS
    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("Fn:" + str(fn) + "[S] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,0],8), end='')
        print()
        print("         ", end="")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,1],8), end='')
        print()
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        print("Fn:" + str(fn) + "[A] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['a_params'][i],8), end='')
        print()
    display.print_line()
  
  def finish_3(results=None):
# HEADER
    display.header_3()
    
# PRINT RSS
    print("# BEST RSS:                  ",g.pfdata['rss']['best'],"") 
    display.print_line()
    
# Best BP
    bp_best = g.pfdata['bp_best']
    
# PRINT VALUES
    print("#        ID   a0       e0       B0       C11      C12      C44 ")
    try:
      for bp_id in bp_best['input'].keys():
        print(  display.pad_r("# Best: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(bp_best['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    try:
      for bp_id in g.bp_known['input'].keys():      
        print(  display.pad_r("# Known: ", 9)
              + display.pad_r(bp_id,4) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['a0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['e0'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['b0_gpa'],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,0],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][0,1],8) + " "
              + display.pad_r(g.bp_known['bp_calculations'][bp_id]['ec_gpa'][3,3],8) + " ")
    except:
      pass
    display.print_line()
    
# PRINT PARAMETERS
    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        print("Fn:" + str(fn) + "[S] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,0],8), end='')
        print()
        print("         ", end="")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['s_nodes'][i,1],8), end='')
        print()
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        print("Fn:" + str(fn) + "[A] ", end='')
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          print("  [" + str(i) + "]", display.pad_r(g.pot_functions['functions'][fn]['a_params'][i],8), end='')
        print()
    display.print_line()
  
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
    print("# Pop Size:    " + str(g.fit['pop_size']) + "      Fresh Size:    " +  str(g.fit['fresh_size']))
    print("# Atoms/sec:        ", end="")
    print(display.pad_l(g.benchmark['atomspersec'], 14), end="")
    print("   # Interactions/sec: ", end="")
    print(display.pad_l(g.benchmark['interationspersec'], 14), end="")
    print("   # Configs/sec:   ", end="")
    print(display.pad_l(g.benchmark['configspersec'], 14), end="")
    print("   # Configs:       ", end="")
    print(display.pad_l(g.benchmark['configs'], 14), end="")
    print()    
    print(display.pad_r_always("# Next Extinction:    " + e_print, 60))    
    display.print_line()
    g.benchmark['configspersec']
    line = ['','','','','']
# Col 1
    line[0] = display.pad_r_always("# Timer:              " + display.pad_l(time.time() - g.times['start'], 16), 60)
    line[1] = display.pad_r_always("# Stage:              " + g.pfdata['stage'], 60)
    line[2] = display.pad_r_always("# RSS Counter:        " + display.pad_l(g.pfdata['rss']['counter'], 16), 60)
    line[3] = display.pad_r_always("# Since Improvement:  " + display.pad_l(g.pfdata['rss']['since_improvement'], 16), 60)
    line[4] = display.pad_r_always("" + e_print, 60)
    
# Col 2
    line[0] = line[0] + "# " + display.pad_r_always("TOP 10", 21)
    
    for n in range(10):
      ln = (n%4) + 1
      if(len(g.top_parameters) > n):
        line[ln] = line[ln] + display.pad_r_always(g.top_parameters[n][0], 20)
      
#if(g.pfdata['top']['filled']):
#  line[ln] = line[ln] + display.pad_r_always(g.top_parameters[n][0], 20)
#else:
#  tn = n + (g.pfdata['top']['size'] - g.pfdata['top']['counter'])
#  if(tn < g.pfdata['top']['size']):
#    line[ln] = line[ln] + display.pad_r_always(g.pfdata['top']['rss'][tn], 20)
    
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

###########################################
#  CLASS pf_dat
###########################################
class pf_data:
  
  def make():
    main.log("Make data structure")

    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    width = potential.parameter_count()
    
    main.log("Pop size:    " + str(pop_size))
    main.log("Fresh size:  " + str(fresh_size))
    main.log("Param width: " + str(width))
    
    top_size = 10
    if(top_size < 10):
      top_size = 10
    top = {'filled': False, 'counter': 0, 'size': top_size, 'rss': numpy.zeros((top_size,),), 'p': numpy.zeros((top_size, width,),),}
    top['rss'][:] = -1.0

    try:
      bp_atom_count = bp.total_atoms
    except:  
      bp_atom_count = 0
      
    try:
      efs_atom_count = efs.total_atoms
    except:  
      efs_atom_count = 0      
      
    d = {
        'stage': 'Not Set', 
        'pool': {'params': None,
                 'pn': 0,
                 'sane_a': 0,
                 'sane_b': 0,
                 'pmax': 10000,
        },
        'params': {'count': width, 
                   'start': numpy.zeros((width,),),
                   'var': numpy.zeros((2,width,),),
                   'best': numpy.zeros((width,),),
                   'pop': numpy.zeros((pop_size, width,),),
                   'pop_rss': numpy.zeros((pop_size,),),
                   'fresh': numpy.zeros((fresh_size, width,),),
                   'fresh_rss': numpy.zeros((fresh_size,),),
                   'children': numpy.zeros((pop_size + 2 * fresh_size, width,),),
                   'children_rss': numpy.zeros((pop_size + 2 * fresh_size,),),
                   'current': numpy.zeros((width,),),
                   'pop_merged': numpy.zeros((2 * pop_size + 3 * fresh_size, width,),),
                   'pop_merged_rss': numpy.zeros((2 * pop_size + 3 * fresh_size,),),
        },
        'rss': {'start': None,
                'current': None,
                'best': None,
                'counter': 0,
                'counter_successful': 0,
                'since_improvement': 0,
        },
        'max_density': {'current': None,
                        'best': None,
        },
        'cycle': {'counter': 0, 'total_cycles': g.fit['cycles'], 'spline_counter': 0,},
        'generation': {'counter': 0, 'total_generations': g.fit['gens'], 'spline_counter': 0,},
        'top': top,
        'extinction': {'event_count': 0, 'every': g.fit['exct_every'], 'counter': 0, 'top_bias': g.fit['exct_top_bias'],},
        'enhance': {'event_count': 0, 'every': g.fit['enhance_every'], 'top': g.fit['enhance_top'], 'counter': 0,},
        'best_hash': None,
        'bp': {},
        'efs': {},
        'bp_known': g.bp_known,
        'efs_known': g.efs_known,
        'bp_best': None,
        'time': {'start': g.times['start'], 'end': None},
    }
    
# SAVE Starting Parameters
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        d['params']['start'][a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        d['params']['start'][a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b    
    
# SAVE Variation
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        d['params']['var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        d['params']['var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        d['params']['var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        d['params']['var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  
    
    return d

###########################################
#  CLASS pf_cycl
###########################################
class pf_cycle:

  def run():  
  
# Increment Cycle Counter
    g.pfdata['cycle']['counter'] += 1
    
# Load from globals
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    
#############################
# Initialise first half
# Search within ranges
#############################
    
    g.pfdata['stage'] = 'Initialising Population'
    
# Add top parameters (if there are any)
    
    cycle_start_parameters = []
    new_p = copy.deepcopy(g.pfdata['params']['start'][:])
    cycle_start_parameters.append(new_p)
    
    for i in range(len(g.top_parameters)):
      new_p = copy.deepcopy(g.top_parameters[i][1][:])
      cycle_start_parameters.append(new_p)

    n = 0
    for pn in range(pop_size):
      loop = True
      while(loop): 
        n = n + 1
        if(n <= len(cycle_start_parameters)):
          g.pfdata['params']['pop'][pn, :] = cycle_start_parameters[n-1][:]
        else:
          g.pfdata['params']['pop'][pn, :] = pf_parameters.random_p()
        
# Try - if it fails or rss == None, try next
        pf_potential.update(g.pfdata['params']['pop'][pn, :])
        rss = pf.get_rss()
        try:
# Update
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['pop_rss'][pn] = rss
        except:
          pass

###################################
# LOOP THROUGH GENERATIONS
###################################
  
    for gen in range(g.fit['gens']):
      g.pfdata['generation']['counter'] += 1
      g.pfdata['stage'] = 'Loop Through Generations - Gen ' + str(g.pfdata['generation']['counter']) 
      
# Run through a generation
      pf_generation.run()
      
  def random_p(c=0.0, m=1.0):
#
    p_count = g.pfdata['params']['count']
    lower = g.pfdata['params']['var'][0,:]
    upper = g.pfdata['params']['var'][1,:]
    range = upper - lower
    
# If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      c = lower + 0.5 * range 
    
# Multiply range
    m_range = m * range
    
# Get random parameters 0 to 1
    r = numpy.random.rand(p_count)
    
# New parameters
    e = 0
    if(g.pfdata['generation']['counter']>0):
      e = g.pfdata['generation']['counter'] - 1
    
    p_new = c + (r - 0.5) * m_range * (g.fit['gen_var_factor'])**(e)
    
# Return
    return p_new
    
  def run_old():  
  
# Increment Cycle Counter
    g.pfdata['cycle']['counter'] += 1
    
# Load from globals
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    
#############################
# Initialise first half
# Search within ranges
#############################
    
    g.pfdata['stage'] = 'Initialising Population - First Half'
    for p in range(pop_size // 2):
      if(p == 0):
        g.pfdata['params']['pop'][p, :] = g.pfdata['params']['start'][:]
        pf_potential.update(g.pfdata['params']['pop'][p, :])
        rss = pf.get_rss()
        g.pfdata['params']['pop_rss'][p] = rss
      else:
        loop = True
        while(loop): 
          g.pfdata['params']['pop'][p, :] = pf_cycle.random_p(0.0, 1.0)
          try:
            pf_potential.update(g.pfdata['params']['pop'][p, :])
            rss = pf.get_rss()
            if(g.pfdata['rss']['current'] is not None):
              loop = False
              g.pfdata['params']['pop_rss'][p] = rss
          except:
            pass

#############################
# Initialise second half
# Search a wider area
#############################
  
    g.pfdata['stage'] = 'Initialising Population - Second Half'
    w = w_start
    w_inc = (w_end - w_start) / (pop_size // 2 - 1)
    for p in range(pop_size // 2, pop_size):
      w = w + w_inc
      loop = True
      while(loop): 
        g.pfdata['params']['pop'][p, :] = pf_cycle.random_p(0.0, w)
        try:
          pf_potential.update(g.pfdata['params']['pop'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['pop_rss'][p] = rss
        except:
          pass
            
###################################
# LOOP THROUGH GENERATIONS
###################################
  
    for gen in range(g.fit['gens']):
      g.pfdata['generation']['counter'] += 1
      g.pfdata['stage'] = 'Loop Through Generations - Gen ' + str(g.pfdata['generation']['counter']) 
      
# Run through a generation
      pf_generation.run()
      
###########################################
#  CLASS pf_potentia
###########################################
class pf_potential:

  def update(p, no_rescale=False):    
  
# Update potential
    a = 0

    for fn in range(len(g.pot_functions['functions'])): 
# Calc b
      b = a + g.pot_functions['functions'][fn]['fit_size']  
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE           
        g.pot_functions['functions'][fn]['s_nodes'][:,1] = p[a:b]
        potential.make_spline_points_inner(fn)       
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
# Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)

# Update a
      a = b     
        
# Rescale embedding data points to cover density range
    potential.rescale_embedding()

# Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials
    
# Save parameters
    potential.parameters = numpy.copy(potential.get_parameters())
    g.pfdata['p']['current'] = numpy.copy(potential.get_parameters())
 
  def take_density(p, p_dens): 
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      b = a + g.pot_functions['functions'][fn]['fit_size'] 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE       
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          p[a:b] = copy.deepcopy(p_dens[a:b])
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC          
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          p[a:b] = copy.deepcopy(p_dens[a:b])
      a = b
    return p

"""
  
# Update potential
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE   
# Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size']  
        
        g.pot_functions['functions'][fn]['s_nodes'][:,1] = p[a:b]
        potential.make_spline_points_inner(fn)
        
# Make Analytic Points
#g.pot_functions['functions'][fn]['s_params'][:] = p[a:b]
# LOAD ORIGINAL
#g.pot_functions['functions'][fn]['points'] = numpy.copy(g.pot_functions['functions'][fn]['points_original'])
# VARY SPLINE
#potential.vary_tabulated_points(fn, p[a:b])
        
# Update a
        a = b   
        
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
# Calc b
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
# Make Analytic Points
        g.pot_functions['functions'][fn]['a_params'][:] = p[a:b]
        potential.make_analytic_points_inner(fn)
        a = b    
"""

###########################################
#  CLASS pf_generatio
###########################################
class pf_generation:

  parents = None
  pop_size = None
  pop_size_half = None
  fresh_size = None

  def run():
    main.log_hr()
    main.log("Generation " + str(g.pfdata['generation']['counter']))
    main.log_hr()

    pf_generation.pop_size = g.fit['pop_size']
    pf_generation.pop_size_half = g.fit['pop_size'] // 2
    pf_generation.fresh_size = g.fit['fresh_size']

# Make array of parents and shuffle
    pf_generation.parents = numpy.arange(pf_generation.pop_size)
    numpy.random.shuffle(pf_generation.parents)
    
    c = 0
# Breed Population
    g.pfdata['stage'] = 'Breed Population'  
    for p in range(pf_generation.pop_size_half):  
      loop = True
      while(loop):
        loop = pf_generation.breed_event(p, c, 'p+p')      
      c = c + 2
       
# Make Fresh
    g.pfdata['stage'] = 'Initialising Fresh Population'   
    pf_generation.make_fresh()
       
# Reshuffle
    numpy.random.shuffle(pf_generation.parents)
    
# Breed Population-Fresh
    g.pfdata['stage'] = 'Breed With Fresh Population'  
    for p in range(pf_generation.fresh_size):  
      loop = True
      while(loop):
        loop = pf_generation.breed_event(p, c, 'p+f')      
      c = c + 2
      
# Merge populations and select top to form next generation
    g.pfdata['stage'] = 'Merge Populations'   
    pf_generation.merge()
    
# Extinction Event
    g.pfdata['stage'] = 'Extinction'   
    pf_extinction.run()
    
# Load best back into population
    g.pfdata['params']['pop'][:, :] = g.pfdata['params']['pop_merged'][0:pf_generation.pop_size, :]
    g.pfdata['params']['pop_rss'][:] = g.pfdata['params']['pop_merged_rss'][0:pf_generation.pop_size]
    
# Enhance
    g.pfdata['stage'] = 'Enhance'   
    pf_enhance.run()

#g.pfdata['generation']['counter']
#pf_setop
    
    cstr = str(g.pfdata['generation']['counter'])
    while(len(cstr)<5):
      cstr = '0' + cstr    
    gen_dir = g.dirs['wd'] + '/fitting/genetic/' +  cstr
    std.make_dir(gen_dir)
  
    potential_output.full(gen_dir)
    
  def breed_event(p, c, opt):
    pc = g.pfdata['params']['count']
    
    pa = pf_generation.parents[p]
    if(opt == 'p+p'):  
      pb = pf_generation.parents[p + pf_generation.pop_size_half]
      pb_array = 'pop'
    if(opt == 'p+f'): 
      pb = p
      pb_array = 'fresh'
    
    ca = c
    cb = c + 1
    
# Breed
    state = pf_generation.get_state()
    for i in range(pc):
      state = pf_generation.get_state(state)
      if(state):
        g.pfdata['params']['children'][ca, i] = g.pfdata['params']['pop'][pa, i]
        g.pfdata['params']['children'][cb, i] = g.pfdata['params'][pb_array][pb, i]
      else:      
        g.pfdata['params']['children'][cb, i] = g.pfdata['params']['pop'][pa, i]
        g.pfdata['params']['children'][ca, i] = g.pfdata['params'][pb_array][pb, i]
      
# Mutate
    g.pfdata['params']['children'][ca, :] = pf_generation.mutate(g.pfdata['params']['children'][ca, :], g.fit['mutate_chance'])
    g.pfdata['params']['children'][cb, :] = pf_generation.mutate(g.pfdata['params']['children'][cb, :], g.fit['mutate_chance'])
    
# No clones - Child A
    for pn in range(pf_generation.pop_size):
      if((g.pfdata['params']['children'][ca, :] == g.pfdata['params']['pop'][pn,:]).all()):
        r = numpy.random.rand(pc)   
        g.pfdata['params']['children'][ca, :] = g.pfdata['params']['children'][ca, :] * g.fit['no_clone_var'] * (0.5 - r)
        break
        
# No clones - Child B
    for pn in range(pf_generation.pop_size):
      if((g.pfdata['params']['children'][cb, :] == g.pfdata['params']['pop'][pn,:]).all()):
        r = numpy.random.rand(pc)   
        g.pfdata['params']['children'][cb, :] = g.pfdata['params']['children'][cb, :] * g.fit['no_clone_var'] * (0.5 - r)
        break

    loop = False 
    
    pf_potential.update(g.pfdata['params']['children'][ca, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
    pf_potential.update(g.pfdata['params']['children'][cb, :])
    rss = pf.get_rss()
    if(rss is None):
      loop = True
      
    if(loop == False):
      g.pfdata['params']['children_rss'][ca] = rss
      g.pfdata['params']['children_rss'][cb] = rss
      
# Loop again or not
    return loop
    
  def get_state(state = None):
    if(state == None):
      if(random.uniform(0.0, 100.0) > 50.0):
        return True
      return False
    else:
      if(random.uniform(0.0, 100.0) > 50.0):
        if(state):
          return False
        return True
      return state
      
  def mutate(params, chance=0.01):
    mutant = pf_cycle.random_p(0.0, g.fit['mutate_scale'])
    for i in range(len(params)):
      if(random.uniform(0.0, 1.0) <= chance):
        params[i] = mutant[i]
    return params   
    
  def make_fresh():    
    w = g.fit['fresh_ws']
    w_inc = (g.fit['fresh_we'] - g.fit['fresh_ws']) / (g.fit['fresh_size'] // 2 - 1)
    for p in range(g.fit['fresh_size']):
      loop = True
      while(loop):
        try:
          rn = min(len(g.top_parameters) - 1, numpy.floor((len(g.top_parameters) + 1) * random.uniform(0.0, 1.0)**3))
          pbest = numpy.copy(g.top_parameters[rn][1])

          g.pfdata['params']['fresh'][p, :] = pf_parameters.random_p(pbest)
          pf_potential.update(g.pfdata['params']['fresh'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['fresh_rss'][p] = rss
        except:
          pass 

# Used to merge parents, fresh and all children into ordered array
  def merge():
    ps = g.fit['pop_size']
    fs = g.fit['fresh_size']
    cs = g.fit['pop_size'] + 2 * g.fit['fresh_size']
    
# Reset array
    g.pfdata['params']['pop_merged'][:,:] = 0.0
    g.pfdata['params']['pop_merged_rss'][:] = 0.0
    g.pfdata['params']['pop_merged'][0:ps,:] = g.pfdata['params']['pop'][:,:]   
    g.pfdata['params']['pop_merged_rss'][0:ps] = g.pfdata['params']['pop_rss'][:] 
    g.pfdata['params']['pop_merged'][ps:ps+fs,:] = g.pfdata['params']['fresh'][:,:]   
    g.pfdata['params']['pop_merged_rss'][ps:ps+fs] = g.pfdata['params']['fresh_rss'][:] 
    g.pfdata['params']['pop_merged'][ps+fs:ps+fs+cs,:] = g.pfdata['params']['children'][:,:]  
    g.pfdata['params']['pop_merged_rss'][ps+fs:ps+fs+cs] = g.pfdata['params']['children_rss'][:] 
        
# Sort
    pf_generation.sort_merged()
    
  def sort_merged():
    sort.sort_1d_dp_asc(g.pfdata['params']['pop_merged_rss'][:])
    g.pfdata['params']['pop_merged_rss'] = sort.apply_keytable_1d_dp(g.pfdata['params']['pop_merged_rss'])
    g.pfdata['params']['pop_merged'] = sort.apply_keytable_2d_dp(g.pfdata['params']['pop_merged'])

###############################################################################################
#
# Delete some day
#
###############################################################################################

  def make_fresh_old():    
    w = g.fit['fresh_ws']
    w_inc = (g.fit['fresh_we'] - g.fit['fresh_ws']) / (g.fit['fresh_size'] // 2 - 1)
    for p in range(g.fit['fresh_size']):
      loop = True
      while(loop):   
        if(p % 2 == 0):
          g.pfdata['params']['fresh'][p, :] = pf_cycle.random_p(g.pfdata['top']['p'][0,:], w)       
        else:  
          r = random.randint(0, 9)
          g.pfdata['params']['fresh'][p, :] = pf_cycle.random_p(g.pfdata['top']['p'][r,:], w)              
        loop = False
        try:
          pf_potential.update(g.pfdata['params']['fresh'][p, :])
          rss = pf.get_rss()
          if(g.pfdata['rss']['current'] is not None):
            loop = False
            g.pfdata['params']['fresh_rss'][p] = rss
        except:
          pass 
      if(p % 2 == 1):
        w = w + w_inc
  
###########################################
#  CLASS pf_extinctio
###########################################
class pf_extinction:

  def run():
  
    g.pfdata['extinction']['counter'] += 1
    if(g.pfdata['extinction']['every'] == g.pfdata['extinction']['counter']):    
      g.pfdata['extinction']['event_count'] += 1
      g.pfdata['extinction']['counter'] = 0
    else:
      return 0
  
    pf_extinction.best = g.pfdata['params']['pop_merged_rss'][0]
    pf_extinction.worst = g.pfdata['params']['pop_merged_rss'][-1]
    
    top_ten = copy.deepcopy(g.pfdata['params']['pop_merged'][0:10])
        
#print("Extinction")
#print(pf_extinction.best, pf_extinction.worst)
    
# Cause extinction of poorly fitting parameters
    for i in range(len(g.pfdata['params']['pop_merged_rss'][:])):
      chance = pf_extinction.chance(g.pfdata['params']['pop_merged_rss'][i])
      r = numpy.random.uniform(0.0,1.0)
      if(r < chance):       
#print("Change", i, g.pfdata['params']['pop_merged_rss'][i])
        loop = True
        while(loop):
          if(numpy.random.uniform(0.0,1.0)<=g.fit['exct_top_bias']):
            rn = 0
          else:
            rn = numpy.random.randint(0,9)
          new_p = pf_cycle.random_p(top_ten[rn], g.fit['exct_var'])
          pf_potential.update(new_p)
          rss = pf.get_rss()
          if(rss != None):
#print(rss, g.pfdata['params']['pop_merged_rss'][i])
            if(rss < g.pfdata['params']['pop_merged_rss'][i]):
              g.pfdata['params']['pop_merged'][i,:] = new_p
              g.pfdata['params']['pop_merged_rss'][i] = rss
            loop = False

# Sort
    pf_extinction.sort_merged()
#print(g.pfdata['params']['pop_merged_rss'])

  def chance(rss):
    return g.fit['exct_factor'] * ((rss - pf_extinction.best) / (pf_extinction.worst - pf_extinction.best)) 

  def sort_merged():
    sort.sort_1d_dp_asc(g.pfdata['params']['pop_merged_rss'][:])
    g.pfdata['params']['pop_merged_rss'] = sort.apply_keytable_1d_dp(g.pfdata['params']['pop_merged_rss'])
    g.pfdata['params']['pop_merged'] = sort.apply_keytable_2d_dp(g.pfdata['params']['pop_merged'])

###########################################
#  CLASS pf_enhanc
###########################################
class pf_enhance:

  def run():
  
    if(g.pfdata['enhance']['top'] == 0):
      return 0
  
    g.pfdata['enhance']['counter'] += 1
    if(g.pfdata['enhance']['every'] == g.pfdata['enhance']['counter']):    
      g.pfdata['enhance']['event_count'] += 1
      g.pfdata['enhance']['counter'] = 0
    else:
      return 0
  
    pf_setop.run()

    return 1

    pop_count = g.pfdata['enhance']['top']
    if(pop_count > g.fit['pop_size']):
      pop_count = g.fit['pop_size']
      
#for p in range(pop_count):
#  print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
    
    for p in range(pop_count):
#print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
      p_diff = pf_parameters.p_diff()
#print(p_diff)
#exit()
      params, rss = gd.opt(pf_enhance.gd_rss, g.pfdata['params']['pop'][p,:], p_diff)
      if(gd.rss_out < g.pfdata['params']['pop_rss'][p]):
        g.pfdata['params']['pop_rss'][p] = gd.rss_out
        g.pfdata['params']['pop'][p,:] = numpy.copy(params)
    
        pf_potential.update(g.pfdata['params']['pop'][p,:])
        rss = pf.get_rss(True)
  
#enhance_freq     #enhance_freq
    
  def process(p_list):
    rss_list = []
    p_list_out = []
    rss_list_out = []
    for n in range(len(p_list)):
      p_before = numpy.copy(p_list[n])
      rss_before = pf.get_rss()
      rss_list.append(rss_before)
      
      p_diff = pf_parameters.p_diff()
      p_after, rss_after = gd.opt(pf_enhance.gd_rss, p_before, p_diff)
      
      if(rss_after < rss_before):
        p_list_out.append(p_after)
        rss_list_out.append(rss_after)
      else:
        p_list_out.append(p_before)
        rss_list_out.append(rss_before)
        
    print(rss_list)
    print(rss_list_out)
    
  def gd_rss(params):
    pf_potential.update(params)
    rss = pf.get_rss(False)
    return rss 

###########################################
#  CLASS pf_parameter
###########################################
class pf_parameters:

  def random_p(c=0.0, m=None):
  
# Multiply the range
    if(m == None):
      m = 1.0

# Randomly pick over sized parameters (defined in input file)
    rn = numpy.random.uniform()    
    b = g.fitting['oversized_parameters'][0]
    prb = 1.0
    for i in range(1, len(g.fitting['oversized_parameters']),1):
      prb = prb * g.fitting['oversized_parameters'][i]
      if(rn<prb):
        m = m * b
      else:
        break

# Get the upper and lower range
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower
    
# If there's no center, take midpoint of upper/lower - else center it on the parameters c
    if(type(c) != numpy.ndarray and c == 0.0):
      center = lower + 0.5 * p_range 
    else:
      center = numpy.copy(c)
    
# Multiply range
    m_range = m * p_range
    
# Get random parameters 0 to 1
    r = numpy.random.rand(p_count)
    
# New parameters
    p_new = center + (r - 0.5) * m_range
    
# Return
    return p_new
    
  def random_variation(p_in, variation):
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower

# p_out
    p_out = numpy.copy(p_in)
    p_out = p_out +  variation * p_range * (numpy.random.rand(p_count) - 0.5)

    return p_out

  def mutate(p_in, variation):
    p_count = g.pfdata['psize']
    lower = g.pfdata['params_var'][0,:]
    upper = g.pfdata['params_var'][1,:]
    p_range = upper - lower

# p_out
    p_out = numpy.copy(p_in)

###############################################################################################
#
# Delete some day
#
###############################################################################################

#e = 0
#if(g.pfdata['generation']['counter']>0):
#  e = g.pfdata['generation']['counter'] - 1
#p_new = c + (r - 0.5) * m_range * (g.fit['gen_var_factor'])**(e)
    
  def get_p_from_pool():
    pmax = g.pfdata['pool']['pmax']
    pn = g.pfdata['pool']['pn'] % pmax
    p = copy.deepcopy(g.pfdata['pool']['params'][pn, :])                     # Copy parameters
    g.pfdata['pool']['params'][pn, :] = pf_parameters.random_p(p, 0.01)      # Disturb
    g.pfdata['pool']['pn'] = g.pfdata['pool']['pn'] + 1                      # Increment
    return p                                                                 # Return
  
  def setup_pool():
    width = potential.parameter_count()
    g.pfdata['pool']['sane_a'] = g.fit['sane_seeds_a']
    g.pfdata['pool']['sane_b'] = g.fit['sane_seeds_b']
    g.pfdata['pool']['pmax'] = g.fit['pool_size']
    
    if(g.pfdata['pool']['pmax'] < (g.pfdata['pool']['sane_a'] + g.pfdata['pool']['sane_b'])):
      g.pfdata['pool']['pmax'] = g.pfdata['pool']['sane_a'] + g.pfdata['pool']['sane_b']
    g.pfdata['pool']['params'] = numpy.zeros((g.pfdata['pool']['pmax'],width,),)
  
  def make_pool(p = None):
    print("Make Pool")
    main.log_title("Make Pool")    
    main.log("pool size:     " + str(g.fit['pool_size']))
    main.log("pop size:      " + str(g.fit['pop_size']))
    main.log("fresh size:    " + str(g.fit['fresh_size']))
    
    if(p == None):
      if(g.pfdata['generation']['counter'] == 0):
        c = 0.0
        main.log("c:    " + str(c))
      else:
        c = g.pfdata['params']['best']
        
        str_c = ''
        for i in range(len(c)):
          str_c = str_c + str(c[i]) + ' '
        main.log("c:    " + str(str_c))
        
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    pmax = g.fit['pool_size']
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    sa = g.pfdata['pool']['sane_a']
    sb = g.pfdata['pool']['sane_b']
    width = potential.parameter_count()
#width = potential.parameter_count()
#g.pfdata['pool']['params'] = numpy.zeros((pmax,width,),)
#g.pfdata['pool']['pmax'] = pmax
  
# Any Density
    if(g.fit['rescale_density'] == 0):
      print("Creating Pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(c, w)
        g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
#print(pn)
    elif(g.fit['rescale_density'] == 1 or g.fit['rescale_density'] == 2):
    
      sane_a = numpy.zeros((sa,width,),)
      sane_b = numpy.zeros((sb,width,),)
      
      print("Make Sane A")
      main.log("make sane set A")
      pn = 0
      while(pn < sa):
        sane_a[pn,:] = pf_cycle.random_p(c)
        a = 0
        for fn in range(len(g.pot_functions['functions'])): 
          b = a + g.pot_functions['functions'][fn]['fit_size'] 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            loop = True
            while(loop):
              pf_potential.update(sane_a[pn,:], True)
              rho = pf_parameters.estimate_density(fn)
              if( rho > 0.0 and rho <=1.0):
                loop = False
              else:
                ptemp = pf_cycle.random_p(c) 
                sane_a[pn,a:b] = numpy.copy(ptemp[a:b])
          a = b
        pn = pn + 1
              
      main.log(str(pn))      
      print("Make Sane B")
      main.log("make sane set B")
      pn = 0        
      while(pn < sb):
        loop = True
        while(loop):
          loop = False
          r = numpy.random.rand(width)
          params = (1.0 + 0.1 * (0.5-r)) * sane_a[pn%sa,:]
          pf_potential.update(params, True)
          for fn in range(len(g.pot_functions['functions'])): 
            if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
              rho = pf_parameters.estimate_density(fn)
              if(not(rho > 0.0 and rho <=1.0)): 
                loop = True             
          sane_b[pn,:] = params
        pn = pn + 1
      main.log(str(pn))
        
      print("Make Pool")
      main.log("make pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - (1+sa+sb))
      pn = 0
      while(pn<pmax):
        p = pf_parameters.random_p(c, w)
        if(g.fit['sane_fraction']>numpy.random.rand()):
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(pf_potential.take_density(p, sane_b[pn%sb,:]))  
        else:
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        
#numpy.savetxt('test.out', g.pfdata['pool']['params'], delimiter=',')
#exit()
        
# Shuffle
    print("Shuffle Pool")
    main.log("shuffle pool")
    numpy.random.shuffle(g.pfdata['pool']['params'])

  def estimate_density(fn):  
#print(g.pot_functions['functions'][fn]['points'][:,1])
    r = numpy.zeros((7,),)
    rn = numpy.zeros((7,),)
    r[0] = 7.48332e0
    r[1] = 6.32456e0
    r[2] = 6.92821e0
    r[3] = 4.89898e0
    r[4] = 5.65686e0
    r[5] = 2.82843e0
    r[6] = 4.0e0
    rn[0] = 48
    rn[1] = 24
    rn[2] = 8
    rn[3] = 24
    rn[4] = 12
    rn[5] = 12
    rn[6] = 6    
    rho = 0.0    
    for i in range(7):
      y = interp.search_x(r[i], g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
      rho = rho + rn[i] * y
#print(rho)
    return rho
    
  def p_diff():
    return (g.pfdata['params']['var'][1,:] - g.pfdata['params']['var'][0,:])

###########################################
#  CLASS pf_seto
###########################################
class pf_setop:

  def run():

    for n in range(g.fit['start_random']):
      p = pf_parameters.random_p()
      pf_potential.update(p)
      rss = pf.get_rss()

    step_enhance_steps = 10
    
# Loop through entire top list
    p_enhance = []
    for i in range(len(g.top_100)):
      p_enhance.append(g.top_100[i][1])

# Make new list of top parameters
    new_list = pf_setop.process_list(p_enhance, step_enhance_steps)

# Empty and re-populate top
    g.top_100 = []
    for i in range(len(new_list)):
      pf_potential.update(new_list[i])
      rss = pf.get_rss()

  def process_list(p_list, steps = 100):

    p_new_list = []

    for n in range(len(p_list)):

      g.pfdata['stage'] = 'STEP_' + str(n+1)

      p_current = numpy.copy(p_list[n])
      rss_current = pf.get_rss()
      rss_start = rss_current

      for n in range(steps):
        p_new = numpy.copy(p_current)
        p_new = p_new + pf_setop.random_step(0.0001)

        pf_potential.update(p_new)
        rss = pf.get_rss(False)

        if(rss < rss_current):
          p_current = p_new
          rss_current = rss
      
      p_new_list.append(p_current)

    return p_new_list

  def random_step(m = 0.01):
# Get random parameters 0 to 1
    return m * (0.5 - numpy.random.rand(g.pfdata['params']['count']))

###########################################
#  CLASS pf_g
###########################################
class pf_gd:

  h = 1.0e-5

  def run():
    
# Start
    p = g.top_parameters[0][1]
    rss = pf_gd.rss(p)
    dp = pf_gd.df(p)
    dp = dp / max(dp)

    gamma = 1.0
    pf_gd.line_search(p, dp, gamma, rss)

    exit()

  def rss(p, top_parameters=False):
    pf_potential.update(p)
    return pf.get_rss(top_parameters)

# Central Difference
  def df(p):
    dp = numpy.zeros((len(p),),)
    for i in range(len(p)):
      p_f = numpy.copy(p)
      p_b = numpy.copy(p)

      p_f[i] = p_f[i] + pf_gd.h
      p_b[i] = p_b[i] - pf_gd.h

      rss_f = pf_gd.rss(p_f)
      rss_b = pf_gd.rss(p_b)
      print(i, rss_b, rss_f)

      dp[i] = (pf_gd.rss(p_f) - pf_gd.rss(p_b)) / (2 * pf_gd.h)
    return dp

  def line_search(p, dp, gamma, best_rss):    
    gamma = 1.0
    while(gamma > 1.0e-24):
      p_test = numpy.copy(p)
      p_test = p_test - gamma * dp
      rss = pf_gd.rss(p_test)
      gamma = 0.25 * gamma
      print(gamma, rss, best_rss)
      if(rss < best_rss):
        p_best = p_test
        best_rss = rss
        best_gamma = gamma
        
###########################################
#  CLASS pf_n
###########################################
class pf_ng:

  R = None
  rss = None

  def run():
    g.pfdata['stage'] = 'NEWTON GAUSS'
    g.pfdata['stage_brief'] = 'NG'

    p = g.top_parameters[0][1]
    pf_ng.set(p)
    pf_ng.residual()

    Jw = len(p)
    Jh = len(pf_ng.R)
    R = numpy.copy(g.rss['residual'])
    h = 1.0e-10
    l = 0.1

    J = numpy.zeros((Jh,Jw),)
    for i in range(len(p)):
      p_b = numpy.copy(p)
      p_b[i] = p_b[i] - h
      pf_potential.update(p_b)      
      rss_b = pf.get_rss(False)
      R_b = numpy.copy(g.rss['residual'])

      p_f = numpy.copy(p)
      p_f[i] = p_f[i] + h
      pf_potential.update(p_f)      
      rss_f = pf.get_rss(False)
      R_f = numpy.copy(g.rss['residual'])

      J[:,i] = (R_f - R_b) / (2.0 * h)

    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_diag = numpy.diag(JTJ)
    JTR = numpy.matmul(JT, R) 

    try:
      dp = numpy.linalg.solve(JTJ + l * JTJ_diag, -JTR)
      p = p - dp
    except:
      pass
 
    pf_potential.update(p)
    pf.get_rss(False)

  def residual():
    rss = pf.get_rss(False)
    pf_ng.R = g.rss['residual']
    pf_ng.rss = rss
    
  def jacobian():
    print("J")

  def set(p):
    pf_potential.update(p)

###########################################
#  CLASS pf_s
###########################################
class pf_sa:

  t = None
  count = 0

  def run():
    pf_sa.count = pf_sa.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 

    if(pf.fit['sa_loops_t'] == 0 or pf.fit['sa_loops_i'] == 0):
      return 0

    g.pfdata['stage'] = 'Simulated Annealing ' + str(pf_sa.count)
    g.pfdata['stage_brief'] = 'SA' + str(pf_sa.count)

# Start - use best
    p = g.pfdata['p']['best']
    pf_potential.update(p)
    rss = pf.get_rss(False)
    p_count = len(p)

# Store Best Parameters
    p_best = numpy.copy(p)
    rss_best = rss

    step = pf.fit['sa_step']  
    for tn in range(pf.fit['sa_loops_t']):
      if(pf.fit['sa_loops_t'] == 1):
        pf_sa.t = pf.fit['sa_temp_start']
        f = 1.0
      else:
        pf_sa.t = pf.fit['sa_temp_start'] - tn *  (pf.fit['sa_temp_start'] - pf.fit['sa_temp_end']) / (pf.fit['sa_loops_t'] - 1)
        f = pf.fit['sa_step_factor']**tn
        g.pfdata['stage'] = 'Simulated Annealing ' +str(pf_sa.count) + ' Loop ' + str(tn + 1)
        rss_start_loop = rss_best
      for n in range(pf.fit['sa_loops_i']):
        loop = True
        while(loop):
          p_new = numpy.copy(p)
          p_new = p_new + f * step * (0.5 - numpy.random.rand(p_count))
          pf_potential.update(p_new)
          rss_new = pf.get_rss(False)
          if(rss_new is not None):
            loop = False
        if(rss_new is not None):
          if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss-rss_new) / pf_sa.t)):
            p = numpy.copy(p_new)
            rss = rss_new
            if(rss_new < rss_best):
              p_best = numpy.copy(p_new)
              rss_best = rss_new
      if(rss_best < rss_start_loop):
        pf_potential.update(p_best)
        rss_new = pf.get_rss(True)

    pf.summary_line("SIMULATED ANNEALING " + str(pf_sa.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
    pf_save.top("SIMULATED_ANNEALING_" + str(pf_sa.count))

###########################################
#  CLASS pf_sa
###########################################
class pf_saa:

  t = None
  count = 0

  def run():
    pf_saa.count = pf_saa.count + 1
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 

#if(pf.fit['sa_loops_t'] == 0 or pf.fit['sa_loops_i'] == 0):
#  return 0

# Start - use best
    p = g.pfdata['p']['best']
    pf_potential.update(p)
    rss = pf.get_rss(False)
    p_count = len(p)

# Store Best Parameters
    p_best = numpy.copy(p)
    rss_best = rss
    g.pfdata['stage'] = 'Simulated Annealing Advanced ' + str(pf_saa.count) + ' Line Search'

    f = 100.0
    f_best = f
    p_grad = pf_pgradient.run(p_best)
    for n in range(100):
      p_new = numpy.copy(p_best)
      p_new = p_new - f * p_grad
      pf_potential.update(p_new)
      rss_new = pf.get_rss(False)
      if(rss_new < rss_best):
        rss_best = rss_new 
        f_best = f
      f = 0.5 * f
 
    g.pfdata['stage'] = 'Simulated Annealing Advanced ' + str(pf_saa.count) + ' Anneal'

    for i in range(10):
      for n in range(50):
        p_new = numpy.copy(p_best)
        p_new = p_new - 0.1 * f_best * p_grad * numpy.random.rand(p_count)   

        pf_potential.update(p_new)
        rss_new = pf.get_rss(False)

        if(rss_new < rss_best):
          p_best = numpy.copy(p_new)
          rss_best = rss_new
      p_grad = pf_pgradient.run(p_best)
   
    exit()

    """
    step = pf.fit['sa_step']  
    for tn in range(pf.fit['sa_loops_t']):
      if(pf.fit['sa_loops_t'] == 1):
        pf_sa.t = pf.fit['sa_temp_start']
        f = 1.0
      else:
        pf_sa.t = pf.fit['sa_temp_start'] - tn *  (pf.fit['sa_temp_start'] - pf.fit['sa_temp_end']) / (pf.fit['sa_loops_t'] - 1)
        f = pf.fit['sa_step_factor']**tn
        g.pfdata['stage'] = 'Simulated Annealing Loop ' + str(tn + 1)
        rss_start_loop = rss_best
      for n in range(pf.fit['sa_loops_i']):
        p_new = numpy.copy(p)
        p_new = p_new + f * step * (0.5 - numpy.random.rand(p_count))
        pf_potential.update(p_new)
        rss_new = pf.get_rss(False)

        if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss-rss_new) / pf_sa.t)):
          p = numpy.copy(p_new)
          rss = rss_new
          if(rss_new < rss_best):
            p_best = numpy.copy(p_new)
            rss_best = rss_new
      if(rss_best < rss_start_loop):
        pf_potential.update(p_best)
        rss_new = pf.get_rss(True)
    """

    pf.summary_line("SIMULATED ANNEALING " + str(pf_random.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])

###########################################
#  CLASS pf_save
###########################################
class pf_saved:

  def run():
  
#pf_pgradient.run()
    g.pfdata['stage'] = 'START - LOAD SAVED'

# Load saved parameters
    try:
      params = pf_top.load(g.dirs['wd'] + '/save', 'top_parameters.txt')    
      for p in params:
        pf_potential.update(p)
        rss = pf.get_rss()
    except:
      pass

###########################################
#  CLASS pf_rando
###########################################
class pf_random:

  count = 0
  time_spent = 0.0

  def run():
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 
    pf_random.count = pf_random.count + 1
    g.pfdata['stage'] = 'RANDOM ' + str(pf_random.count)
    g.pfdata['stage_brief'] = 'R' + str(pf_random.count)

# START PARAMETER
# Try randomly generated parameters
    for n in range(pf.fit['random_size']):
      p = pf_parameters.random_p(0.0)
      pf_potential.update(p)
      rss = pf.get_rss()

    pf_random.time_spent = pf_random.time_spent + (time.time() - t_start)

    pf.summary_line("RANDOM SEARCH " + str(pf_random.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
   
    pf_save.top("RANDOM_SEARCH_" + str(pf_random.count))

###########################################
#  CLASS pf_to
###########################################
class pf_top:

  def save(dir, filename):

    print("Saving Top Parameters")
    print(dir + "/" + filename)
    std.make_dir(dir)
    fh = open(dir + "/" + filename, 'w')
    for i in range(len(g.top_parameters)):
      i_str = str(i+1)
      while(len(i_str)<8):
        i_str = "0" + i_str
      fh.write("---START---" + "\n")
      fh.write(i_str + "\n")
      fh.write("PARAMETERS \n")
      fh.write(str(len(g.top_parameters[i][1])) + "\n")
      for j in range(len(g.top_parameters[i][1])):
        fh.write(str(g.top_parameters[i][1][j]) + "\n")
      fh.write("RSS " + str(g.top_parameters[i][0]) + '\n')
      fh.write("BP MaxDensity " + str(g.top_parameters[i][4]) + '\n')
      fh.write("EFS MaxDensity " + str(g.top_parameters[i][5]) + '\n')
      try:
        bp_calculations = g.top_parameters[i][3]['bp_calculations']
        for bi in range(len(bp_calculations)):
          fh.write("a0 " + str(bp_calculations[bi]['a0']) + "\n")
          fh.write("e0 " + str(bp_calculations[bi]['e0']) + "\n")
          fh.write("b0 " + str(bp_calculations[bi]['b0']) + "\n")
      except:
        pass
      fh.write("---END---" + "\n")

    fh.close()

  def load(dir, filename):
    if(int(g.fitting['load_top_parameters']) == 0):
      return []
    print(g.fitting['load_top_parameters'])
    top_p = []
    fh = open(dir + "/" + filename, 'r')
    
    f = []
    for line in fh:
      f.append(line.strip())
    fh.close()

    i = 0
    while(i<len(f)):
      if(f[i] == "---START---"):
        i = i + 3
        pn = int(f[i])
        p = numpy.zeros((pn,),)
        for n in range(pn):
          i = i + 1
          p[n] = float(f[i])
      elif(f[i] == "---END---"):
        top_p.append(p)
      i = i + 1
      if(len(top_p) == g.fitting['load_top_parameters']):
        break

    return top_p 

###########################################
#  CLASS pf_pgradien
###########################################
class pf_pgradient:

  def run(p):
    pf_potential.update(p)
    rss = pf.get_rss(False)
    grad = numpy.zeros((len(p),), dtype=numpy.float64,)

    for i in range(len(p)):
# Forward
      p_f = numpy.copy(p)
      h = 1.0e-10 * p_f[i]
      if(p_f[i] == 0.0):
        h = 1.0e-10
      p_f[i] = p_f[i] + h
      pf_potential.update(p_f)
      rss_f = pf.get_rss(False)
# Backward
      p_b = numpy.copy(p)
      h = 1.0e-10 * p_b[i]
      if(p_b[i] == 0.0):
        h = 1.0e-10
      p_b[i] = p_b[i] - h
      pf_potential.update(p_b)
      rss_b = pf.get_rss(False)
      
      grad[i] = (rss_f - rss_b) / (2.0 * h)

    return grad

###########################################
#  CLASS pf_displa
###########################################
class pf_display:

  display_header_set = False
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
      pf_display_simple.output()
#if(output == 2):
#  display.output_2()
#if(output == 3):
#  display.output_3()
#  #print(g.rss)

  @staticmethod
  def pad_r(inp, p=12, r=None):
    if(inp == None):
      return ""  
    if(r is not None):  
      try:
        inp = numpy.round(inp, r)
      except:
        pass
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]

  @staticmethod
  def pad_l(inp, p=12, r=None):
    if(inp == None):
      return ""  
    if(r is not None):  
      try:
        inp = numpy.round(inp, r)
      except:
        pass  
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
    
###########################################
#  CLASS pf_display_simpl
###########################################
class pf_display_simple:

  def output():
    if(g.pfdata['last_stage'] == None or g.pfdata['last_stage'] != g.pfdata['stage']):
      g.pfdata['last_stage'] = g.pfdata['stage']
      print("#####################################")
      print(g.pfdata['stage'])
      print("#####################################")
    
    if(g.pfdata['rss']['current'] == None):
      rss_current = "Error"
    else:
      rss_current = '{:12.4e}'.format(g.pfdata['rss']['current'])

    pf_time = '{:6.3f}'.format(time.time() - pf.start_time)
    rss_best = '{:12.4e}'.format(g.pfdata['rss']['best'])

    sb = g.pfdata['stage_brief']
    while(len(sb)<8):
      sb = sb + " "

    print(sb, end = " ")
    print(pf_display.pad_l(pf_time, 10), end="    ")
    print(pf_display.pad_r(rss_best, 14), end=" ")
    print(pf_display.pad_r(rss_current, 14), end=" ")
    print()

###########################################
#  CLASS pf_ini
###########################################
class pf_init:

  def run():  

# If a spline fit, convert into a spline with the set number of nodes
    for fn in range(len(g.pot_functions['functions'])):       
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     
        potential.vary_tabulated_points(fn)

# Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    efs_calc.set_weights()             # Set weightings
    potential.efs_add_potentials()     # Load potentials
    configs.efs_add_config()           # Add configs
    
# Setup BP
    bp_calc.init()
    bp_calc.set_weights()
    potential.bp_add_potentials()
    b_props.bp_add()
    bp_calc.get_known()

    g.pfdata = {}
    
    g.pfdata['stage'] = ''
    g.pfdata['last_stage'] = None

# Set up RSS
    g.pfdata['rss'] = {'start': None,
                       'current': None,
                       'best': None,
                       'counter': 0,
                       'counter_successful': 0,
                       'since_improvement': 0,}

    g.pfdata['psize'] = potential.parameter_count()
    g.pfdata['p'] = {'current': None, 'best': None,}

    g.pfdata['bp'] = {'current': None, 'best': None,}

    g.pfdata['max_density'] = {'bp_current': None, 'efs_current': None, 'bp_best': None, 'efs_best': None,}
    
# SAVE Starting Parameters
    g.pfdata['params_start'] = numpy.zeros((g.pfdata['psize'],),)
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        g.pfdata['params_start'][a:b] = numpy.zeros((g.pot_functions['functions'][fn]['fit_size'],),)
        a = b
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        g.pfdata['params_start'][a:b] = g.pot_functions['functions'][fn]['a_params'][:]
        a = b    
    
# SAVE Variation
    g.pfdata['params_var'] = numpy.zeros((2,g.pfdata['psize'],),)
    a = 0
    for fn in range(len(g.pot_functions['functions'])):        
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # TABULATED      
        b = a + g.pot_functions['functions'][fn]['fit_size']     
        g.pfdata['params_var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        g.pfdata['params_var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b 
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC  
        b = a + g.pot_functions['functions'][fn]['fit_size'] 
        g.pfdata['params_var'][0,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][0,:]  # Lower
        g.pfdata['params_var'][1,a:b] = g.pot_functions['functions'][fn]['fit_parameters'][1,:]  # Upper
        a = b  

###########################################
#  CLASS pf_sav
###########################################
class pf_save:
 
  count = 0

  def top(name):
    pf_save.count = pf_save.count + 1
 
    c = str(pf_save.count)
    while(len(c)<4):
      c = '0' + c
    name = c + '_' + name

    dir_out = g.dirs['wd'] + '/fitting/saved_potentials/' + name + '/'
    plot_dir = dir_out + 'plots/'
    pot_dir = dir_out + 'pots/'
    std.make_dir(dir_out)
    std.make_dir(plot_dir)
    std.make_dir(pot_dir)
    
    p = numpy.copy(g.top_parameters[0][1])
    pf_potential.update(p)
    rss = pf.get_rss()

    fh = open(dir_out + 'summary.txt', 'w')
    fh.write("RSS:  " + '\n')
    fh.write(str(g.top_parameters[0][0]) + '\n')
    fh.write("" + '\n')
    fh.write("P:  " + '\n')
    for pn in p:
      fh.write(str(pn) + '\n')
    fh.write("" + '\n')
    fh.write("EFS Density:" + '\n')
    fh.write(str(g.top_parameters[0][4]) + '\n')
    fh.write("" + '\n')
    fh.write("BP Density:" + '\n')
    fh.write(str(g.top_parameters[0][5]) + '\n')
    fh.write("" + '\n')
    fh.write("" + '\n')
    for i in range(len(g.top_parameters[0][3]['bp_calculations'])):
      
      fh.write("Known:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.bp_known['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.bp_known['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.bp_known['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.bp_known['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write('\n')

      fh.write("Calculated:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.top_parameters[0][3]['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write("" + '\n')

    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.write("" + '\n')

    fh.write("RSS Details:  " + '\n')
    fh.write("===============================  " + '\n')
    fh.write("Total:      " + str(g.top_parameters[0][0]) + '\n')
    fh.write("a0:         " + str(bp.rss_by_type[0]) + '\n')
    fh.write("e0:         " + str(bp.rss_by_type[1]) + '\n')
    fh.write("b0:         " + str(bp.rss_by_type[2]) + '\n')
    fh.write("ec:         " + str(bp.rss_by_type[3]) + '\n')
    fh.write("g:          " + str(bp.rss_by_type[4]) + '\n')
    fh.write("e:          " + str(bp.rss_by_type[5]) + '\n')
    fh.write("v:          " + str(bp.rss_by_type[6]) + '\n')
    fh.write("a0 (w):     " + str(bp.rss_by_type_w[0]) + '\n')
    fh.write("e0 (w):     " + str(bp.rss_by_type_w[1]) + '\n')
    fh.write("b0 (w):     " + str(bp.rss_by_type_w[2]) + '\n')
    fh.write("ec (w):     " + str(bp.rss_by_type_w[3]) + '\n')
    fh.write("g (w):      " + str(bp.rss_by_type_w[4]) + '\n')
    fh.write("e (w):      " + str(bp.rss_by_type_w[5]) + '\n')
    fh.write("v (w):      " + str(bp.rss_by_type_w[6]) + '\n')

    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.close()

    potential.plot_python_potentials(plot_dir)
    potential_output.full(pot_dir)

    pf_top.save(dir_out, 'top_parameters.txt')

###########################################
#  CLASS pf_fina
###########################################
class pf_final:

  def run():

    p = g.pfdata['p']['best']
    pf_potential.update(p)

# Run EFS and BP
    rss = rss_calc.run_calc()

    print('') 
    print('') 
    print('')     
    print('###########################################################################') 
    print('###########################################################################') 
    print('FINAL CALCULATION WITH BEST PARAMETERS') 
    print('###########################################################################') 
    print('###########################################################################') 
    print('')     

    b_props.bp_output_terminal()

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
      print('')   
 
      print('BP  ' + str(bp_id)) 
      print('================') 

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
    print('RSS: ' + str(rss))
    print('')
    print('')
    print('Max Density:', bp.max_density, efs.max_density)
    print('')
    print('')
#print(g.rss)
    print('')
    
###########################################
#  CLASS tria
###########################################
class trial:

  def run():  
  
    print("Trial") 

    print("Saving Configs")
    for cn in range(len(g.configs['configs'])):
      configs.save(cn)
      
    print("Save Potentials")
    potential.save_potential()
      
    potential.plot_python_potentials()

###########################################
#  CLASS relax_cal
###########################################
class relax_calc:

  def run():  
  
    print("Relax") 
    
    relax.init()
    potential.relax_add_potentials()
    
    c = g.configs['configs'][0]  
#relax_calc.randomise_coords(c)
    relax.add_config( 
                     7.0,
                     c['alat'], 
                     c['uv'], 
                     c['coords_label_id'], 
                     c['coords']
                    )  
               
    relax.run(25, 0.07, -0.9)  
               
    relax_calc.save_xyz()
    
  def md():  
    relax.run_md(1000, 0.01)                
    relax_calc.save_xyz_md()
    
  def randomise_coords(c):     
    for ci in range(len(c['coords'])):
      r = numpy.random.rand(3)
      c['coords'][ci, :] = c['coords'][ci, :]  + (0.5- r) * 0.04
    
  def save_xyz():
    dir = g.dirs['wd'] + "/relax"
    std.make_dir(dir)
    fh = open(dir + '/relaxed.xyz', 'w')
    fh.write(str(relax.atom_count) + "\n")
    fh.write("\n")
    for i in range(relax.atom_count):
      fh.write(str(labels.get(relax.labels[i])) + " ")
      fh.write(str(relax.coords[i, 3]) + " ")
      fh.write(str(relax.coords[i, 4]) + " ")
      fh.write(str(relax.coords[i, 5]) + " ")
      fh.write("\n")      
    fh.close() 
    fh = open(dir + '/relaxed_crystal.xyz', 'w')
    fh.write(str(relax.atom_count) + "\n")
    fh.write("\n")
    for i in range(relax.atom_count):
      fh.write(str(labels.get(relax.labels[i])) + " ")
      fh.write(str(relax.coords[i, 0]) + " ")
      fh.write(str(relax.coords[i, 1]) + " ")
      fh.write(str(relax.coords[i, 2]) + " ")
      fh.write("\n")      
    fh.close() 
  
  def save_xyz_md():
    dir = g.dirs['wd'] + "/md"
    std.make_dir(dir)
    
    fh = open(dir + '/md.xyz', 'w')
    for n in range(relax.md_steps + 1):
      fh.write(str(relax.atom_count) + "\n")
      fh.write(str(n) + "\n")
      for i in range(relax.atom_count):
        fh.write(str(labels.get(relax.labels[i])) + " ")
        fh.write(str(relax.md_xyz[n, i, 0]) + " ")
        fh.write(str(relax.md_xyz[n, i, 1]) + " ")
        fh.write(str(relax.md_xyz[n, i, 2]) + " ")
        fh.write("\n")
      
    fh.close()           
               
#relax.make_nl()
#relax.force_calc()
#print(relax.nl_count)
#print(relax.config_forces[0:relax.atom_count, :])
    
# Setup ES
#es.init()
#potential.es_add_potentials()
#es_calc.es_add()

###########################################
#  CLASS main(
###########################################
class main():

  log_info_count = 0

  def start():
  
# RECORD START TIME
    g.times['start'] = time.time()
    
    now = datetime.datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print(date_time)	
  
# OPEN LOG
    main.log_start()

    if(len(sys.argv)>1):
      run_program = False
      try:
        g.inp = read_config.read_file(sys.argv[1])
        main.log('Loaded: ' + str(sys.argv[1]))
        run_program = True
      
# Copy input
#std.copy(sys.argv[0], g.dirs['input'])
#std.copy(sys.argv[1], g.dirs['input'])
      except:
        main.log('Unable to load, exiting\n')
      
# RUN
      if(run_program):
        eampa.run()
        
# End and close log
      main.end()
      
  def log_start():
    fh = open('job.log', 'w')
    fh.close()
    main.log_hr()
    main.log(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
    main.log_hr()
    main.log_br()
    main.log_info()
    main.log('Script: ' + str(sys.argv[0]))
    
  def log_end():
    main.log_br()
    main.log_hr()
    main.log('Duration: ' + str(g.times['end'] - g.times['start']))
    main.log_hr() 
      
  def log(line='', end='\n'):
    fh = open('job.log', 'a') 
    main.log_info_count = main.log_info_count + 1
    if(main.log_info_count == 100):
      main.log_info_count = 0
      main.log_info()  
    lines=line.split("\n")
    for line in lines:
      fh.write(str("{:.3E}".format(time.time() - g.times['start'])) + ' ###   ' + line + end) 
    fh.close()
      
  def log_hr(): 
    fh = open('job.log', 'a') 
    fh.write('#############################################################################\n')
    fh.close() 
    
  def log_br():
    fh = open('job.log', 'a')  
    fh.write('\n') 
    fh.close()
        
  def log_info():
    fh = open('job.log', 'a')  
    fh.write('####TIME#############LOG#####################################################\n') 
    fh.close()
    
  def log_title(title=''): 
    fh = open('job.log', 'a') 
    fh.write('#############################################\n') 
    titles = title.split("\n")
    for title in titles:
      fh.write('  ' + title + '\n') 
    fh.write('  Time: ' + str("{:.3E}".format(time.time() - g.times['start'])) + '\n') 
    fh.write('#############################################\n')
    fh.close() 
  
  def end(): 
# CLOSE LOG
    g.times['end'] = time.time()
    main.log_end()
    exit()

# Run
main.start()

###########################################
###########################################
#  MAIN
###########################################
###########################################

class main():



  log_info_count = 0



  def start():

  

# RECORD START TIME
    g.times['start'] = time.time()

    

    now = datetime.datetime.now()

    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    print(date_time)	

  

# OPEN LOG
    main.log_start()



    if(len(sys.argv)>1):

      run_program = False

      try:

        g.inp = read_config.read_file(sys.argv[1])

        main.log('Loaded: ' + str(sys.argv[1]))

        run_program = True

      

# Copy input
#std.copy(sys.argv[0], g.dirs['input'])
#std.copy(sys.argv[1], g.dirs['input'])
      except:

        main.log('Unable to load, exiting\n')

      

# RUN
      if(run_program):

        eampa.run()

        

# End and close log
      main.end()

      

  def log_start():

    fh = open('job.log', 'w')

    fh.close()

    main.log_hr()

    main.log(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))

    main.log_hr()

    main.log_br()

    main.log_info()

    main.log('Script: ' + str(sys.argv[0]))

    

  def log_end():

    main.log_br()

    main.log_hr()

    main.log('Duration: ' + str(g.times['end'] - g.times['start']))

    main.log_hr() 

      

  def log(line='', end='\n'):

    fh = open('job.log', 'a') 

    main.log_info_count = main.log_info_count + 1

    if(main.log_info_count == 100):

      main.log_info_count = 0

      main.log_info()  

    lines=line.split("\n")

    for line in lines:

      fh.write(str("{:.3E}".format(time.time() - g.times['start'])) + ' ###   ' + line + end) 

    fh.close()

      

  def log_hr(): 

    fh = open('job.log', 'a') 

    fh.write('#############################################################################\n')

    fh.close() 

    

  def log_br():

    fh = open('job.log', 'a')  

    fh.write('\n') 

    fh.close()

        

  def log_info():

    fh = open('job.log', 'a')  

    fh.write('####TIME#############LOG#####################################################\n') 

    fh.close()

    

  def log_title(title=''): 

    fh = open('job.log', 'a') 

    fh.write('#############################################\n') 

    titles = title.split("\n")

    for title in titles:

      fh.write('  ' + title + '\n') 

    fh.write('  Time: ' + str("{:.3E}".format(time.time() - g.times['start'])) + '\n') 

    fh.write('#############################################\n')

    fh.close() 

  

  def end(): 

# CLOSE LOG
    g.times['end'] = time.time()

    main.log_end()

    exit()



# Run
main.start()
