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
from f2py_lib.f_fnc import fnc
from f2py_lib.f_interp import interp
import copy
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
from f2py_lib.f_spline import spline
from f2py_lib.f_bp import polyfit
from f2py_lib.f_es import es
from f2py_lib.f_sorting import sort
import random
import hashlib
from f2py_lib.f_relax import relax

###########################################
#  CLASS
###########################################
class g: 
  
  dirs = {
         'wd': 'wd',
         }
  
  sub_dirs = {
         'log': 'log',  
         'output': 'output',   
         'results': 'results',   
         'plots': 'plots',  
         'eos': 'plots/eos', 
         'ec': 'plots/ec', 
         'pots': 'plots/pots',  
         'fitting': 'fitting',  
         'configs': 'configs',  
         'input': 'input', 
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
    
# Copy Input
    eampa.copy_input()   
    
# Read Input
    read_input.run()    
    
# Set memory
    memory.run() 
    
# Mask
    mask.process()
    
# Load potentials
    potential.load()
    
# Load configs
    configs.load()
    
# Energy Adjustments (dft)
    e_adjust.load()
    
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
  
###########################################
#  CLASS mas
###########################################
class mask:

  def process():    
    m = None
    for k in g.inp.keys():
      if(k.upper() == "MASK"):
        m = k    
    
    if(m == None):
      return ''
        
    g.mask = {}
    for k in g.inp[m].keys():
      g.mask[k.upper()] = g.inp[m][k].upper()
      
  def get(label):
    label = label.upper()
    if(label in g.mask.keys()):
      return g.mask[label]
    return label
    
###########################################
#  CLASS label
###########################################
class labels:

  @staticmethod
  def add(label):
    label = label.upper()
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
    
###########################################
#  CLASS read_inpu
###########################################
class read_input:
 
  def run():
    read_input.run_type()
    read_input.rss_weights()
    read_input.fit()
    read_input.fit_results()

# READ TYPE
  def run_type():
  
# DEFAULT
    g.run_type = 'efs'
    
# TRY READING
    try:
      g.run_type = g.inp['run']['type'].lower().strip()
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
    }
    
# TRY READING
    for k in g.rss_weights.keys():
      try:
        g.rss_weights[k] = float(g.inp['rss_weights'][k])
      except:
        pass
# READ
      
  def fit():
  
# DEFAULT
    g.fit = {
    'cycles': 1,
    'gens': 4,
    'spline_cycles': 0,
    'spline_gens': 0,
    'pop_size': 20,
    'fresh_size': 10,
    'exct_factor': 0.5,
    'exct_every': 5,
    'exct_var': 0.1,
    'exct_top_bias': 0.5,
    'rescale_density': 0,
    'rescale_min': 0.3,
    'rescale_max': 0.9,
    'rescale_default': 0.6,
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
    }

# TRY READING
    for k in g.fit.keys():
      try:
        g.fit[k] = g.inp['fit'][k]
      except:
        pass
        
# POP SIZE - must be even
    if(g.fit['pop_size'] < 2):
      g.fit['pop_size'] = 2
    if(g.fit['pop_size'] % 2 != 0):
      g.fit['pop_size'] = g.fit['pop_size'] + 1
      
# FRESH SIZE - must be even
    if(g.fit['fresh_size'] < 2):
      g.fit['fresh_size'] = 2
    if(g.fit['fresh_size'] % 2 != 0):
      g.fit['fresh_size'] = g.fit['fresh_size'] + 1
      
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
      
###########################################
#  CLASS memor
###########################################
class memory:

  def run():
  
    g.memory = {}

# Defaults
    g.memory['bp'] = {}    
    g.memory['efs'] = {}
      
# For BP module
    try:
      mem = str(g.inp['mem']['bp'])
    except:
      mem = "500MB"
    g.log_fh.write(mem + '\n')
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
    g.log_fh.write(mem + '\n')
    mem = std.mem_value(mem)  
    g.memory['efs']['mem'] = mem    
    g.memory['efs']['c'] = int(1 * (mem / 7840))
    g.memory['efs']['g'] = int(12 * (mem / 7840))
    g.memory['efs']['nl'] = int(100 * (mem / 7840))

#print(g.memory)

###########################################
#  CLASS potentia
###########################################
class potential:

  def run():

    print("Potential") 
  
    efs.init()                           # Initialise (allocate arrays)
    potential.efs_add_potentials()       # Load potentials
    potential.plot_fortran_potentials()
    potential.plot_python_potentials()
    
  def load():
    if(g.inp['potential']['dir'].strip() == ""):
      g.inp['potential']['pot_file'] = g.inp['potential']['index_file']
    else:
      g.inp['potential']['pot_file'] = g.inp['potential']['dir'] + "/" + g.inp['potential']['index_file']
    pot_file = g.inp['potential']['pot_file']
    if(not os.path.isfile(pot_file)):
      return False
# Read potential index
    potential.read_potential(pot_file)
    potential.load_tabulated()
    potential.make_tabulated_points()
    potential.load_analytic()
    potential.make_analytic_points()
    potential.load_fit_data()
    potential.pf_output()
    potential.make_copies()
       
    return True
  
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
    'a_type': '',
    'f': None,
    'a_params': None,
    'a_params_fixed': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
    'points_original': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
    'fit_file': None,
    'fit_type': None,         # 1 spline, 2 analytic
    'fit_parameters': None,
    'fit_size': None,
    'fit_mult': None,
    }

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
        
    index = std.config_file_to_list(file_name)  
    pot = potential.pot_function()
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
        pot['a_text'] = label_str
        pot['a'] = label_id
        if(len(row) > 2):
          label_str, label_id = labels.add(row[2])
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
      elif(len(row) > 1 and row[0].upper() == "ZOOR"):
        val = row[1].upper() 
        pot['zoor'] = 0  
        if(val[0] == "T" or val[0] == "1" or val[0] == "Y"):
          pot['zoor'] = 1  
      elif(len(row) > 0 and row[0].upper() == "END"):
        if(pot['a_text'] != None and pot['f_group_label'] != None):
          label, id = labels.add_group(pot['a_text'], pot['f_group_label'], len(g.pot_functions['functions']))
          pot['f_group'] = id
      
        g.pot_functions['functions'].append(pot)
        pot = potential.pot_function()
    
# READ ZBL DATA
    
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
    
  @staticmethod
  def load_fit_data(): 
# READ FIT DATA
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['fit_file'] != None):
        params = None
        fit_file = g.pot_functions['pot_dir'] + '/' + g.pot_functions['functions'][fn]['fit_file']
        if(os.path.isfile(fit_file)):
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
            if(line[0:2] == "#M"):
              g.pot_functions['functions'][fn]['fit_mult'] = numpy.zeros((2,),)
              """
              try:
                g.pot_functions['functions'][fn]['fit_mult'][0] = float(f[1])
                g.pot_functions['functions'][fn]['fit_mult'][1] = float(f[2])
              except:  
                g.pot_functions['functions'][fn]['fit_mult'][0] = 0.1
                g.pot_functions['functions'][fn]['fit_mult'][1] = 10.0
              """
                
          fh.close()
        if(g.pot_functions['functions'][fn]['fit_type'] != None and type(params_lower) == list):
        
# Analytic fitting
          if(g.pot_functions['functions'][fn]['fit_type'] == 2 and
             g.pot_functions['functions'][fn]['function_type'] == 2):
            g.pot_functions['functions'][fn]['fit_size'] = len(g.pot_functions['functions'][fn]['a_params'])
            g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((2,len(g.pot_functions['functions'][fn]['a_params']),),)     
            
            for i in range(g.pot_functions['functions'][fn]['fit_size']):
              g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
              g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])
              
# Spline fitting
          else:
            s = len(params_lower)
            g.pot_functions['functions'][fn]['fit_type'] = 1
            g.pot_functions['functions'][fn]['fit_size'] = s            
            g.pot_functions['functions'][fn]['fit_parameters'] = numpy.zeros((2,s,),) 
            
            for i in range(len(params_lower)):
              g.pot_functions['functions'][fn]['fit_parameters'][0,i] = float(params_lower[i])
              g.pot_functions['functions'][fn]['fit_parameters'][1,i] = float(params_upper[i])

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
      if(row.strip().upper() == '#A'):
        return 'A'
      return 'T'
      
  @staticmethod
  def pf_output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'pot.dat', 'w')
      
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
    if(dir == None):
      dir = g.dirs['pots']
      
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
      
  def plot_fortran_potentials(dir = None):  
    
    if(dir == None):
      dir = g.dirs['pots']
    
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
# F2PY functions
###########################################################

  def efs_add_potentials():
  
# Clear
    efs.clear_potentials()
    
# ADD POTENTIALS
    for pn in range(len(g.pot_functions['functions'])):
      pf = g.pot_functions['functions'][pn]
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
        
# ADD ZBL
    for zbl in g.pot_functions['zbl']:
      efs.add_zbl(
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
    efs.set_potentials()   
    
  def es_add_potentials():
  
# Clear
    es.clear_potentials()
    
# ADD POTENTIALS
    for pf in g.pot_functions['functions']:
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
    for pf in g.pot_functions['functions']:
      if(pf['f_on']):
        if(pf['f_type_id'] == 1):
          bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['b'],  
                            pf['r_cut'], 
                            pf['points'] 
                           )
        elif(pf['f_type_id'] > 1):
          bp.add_potential(
                            pf['f_type_id'],
                            pf['a'], 
                            pf['f_group'],  
                            pf['r_cut'], 
                            pf['points'] 
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
#  CLASS potential_outpu
###########################################
class potential_output:

  @staticmethod
  def full():
  
# Make output directory
    g.dirs['results'] = g.dirs['wd'] + '/' + g.fit_results['results_dir'] 
    std.make_dir(g.dirs['results'])
    
    potential_output.best_parameters() 
    potential_output.eampa()
    potential_output.data_file()
    potential_output.dl_poly()
    potential_output.plots()
    
  def best_parameters():
    dir_out = g.dirs['results'] + '/best_parameters'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/best_parameters.txt', 'w')    
    for fn in range(len(g.pot_functions['functions'])):          
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # SPLINE
        fh.write("Fn: " + str(fn) + "   Type: spline")
#for i in range(g.pot_functions['functions'][fn]['fit_size']):
#  fh.write("[" + str(i) + "]" + display.pad_r(g.pot_functions['functions'][fn]['fit_parameters'][0,i],8))
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC
        fh.write("Fn:" + str(fn) + "[A] ")
        for i in range(g.pot_functions['functions'][fn]['fit_size']):
          fh.write("  [" + str(i) + "] " + display.pad_r(g.pot_functions['functions'][fn]['a_params'][i],8) + " ")
        fh.write("\n")
    fh.close()
       
  @staticmethod
  def eampa():
    dir_out = g.dirs['results'] + '/pot_save'
    std.make_dir(dir_out)
    
    fh = open(dir_out + '/out.pot', 'w')
    fh.write('POTNAME ' + '\n')  
    fh.write('\n')  
    fh.write('\n')  
    
    for fn in range(len(g.pot_functions['functions'])): 
      f = g.pot_functions['functions'][fn]
#print("F ", fn)
#print(f)
#print()
#print()
# File Name
      pf_name = f['f_type'] + '_' + f['a_text'].strip()
      if(f['b_text'].strip() != ''):
        pf_name += '_' + f['b_text'].strip()
      if(str(f['f_group']).strip() != ''):
        pf_name += '_' + str(f['f_group'])
      pf_name = pf_name.lower() + '.pot'
      
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
    fh.close()
  
  @staticmethod
  def data_file():
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
  def dl_poly():
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
  def plots(): 
    dir_out = g.dirs['results'] + '/plots/python'
    std.make_dir(dir_out)   
    potential.plot_python_potentials(dir_out)
    dir_out = g.dirs['results'] + '/plots/fortran'
    std.make_dir(dir_out)   
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
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        min_rho = min(g.pot_functions['functions'][fn]['points'][:,1])
        g.pot_functions['functions'][fn]['points'][:,1] = g.pot_functions['functions'][fn]['points'][:,1] - min_rho 
        rho = rescale_density.estimate_density(fn) 
        if(rho < g.fit['rescale_min'] or rho > g.fit['rescale_max']):
          g.pot_functions['functions'][fn]['points'][:,1:3] = (g.fit['rescale_default'] / rho) * g.pot_functions['functions'][fn]['points'][:,1:3]
          rho = rescale_density.estimate_density(fn)

  def estimate_density(fn):
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
    return rho
    
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
      print(file) #DEL
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
# Check for mask
          l = mask.get(f[0])        
        
          label_str, label_id = labels.add(l)
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
        fh.write('uv_prim   ' + str(c['uv_prim'][0,0]) + ' ' + str(c['uv_prim'][0,1]) + ' ' + str(c['uv_prim'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][1,0]) + ' ' + str(c['uv_prim'][1,1]) + ' ' + str(c['uv_prim'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv_prim'][2,0]) + ' ' + str(c['uv_prim'][2,1]) + ' ' + str(c['uv_prim'][2,2]) + ' ' + ' \n')
        fh.write('uv        ' + str(c['uv'][0,0]) + ' ' + str(c['uv'][0,1]) + ' ' + str(c['uv'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][1,0]) + ' ' + str(c['uv'][1,1]) + ' ' + str(c['uv'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['uv'][2,0]) + ' ' + str(c['uv'][2,1]) + ' ' + str(c['uv'][2,2]) + ' ' + ' \n')
        fh.write('stress    ' + str(c['stress'][0,0]) + ' ' + str(c['stress'][0,1]) + ' ' + str(c['stress'][0,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][1,0]) + ' ' + str(c['stress'][1,1]) + ' ' + str(c['stress'][1,2]) + ' ' + ' \n')
        fh.write('          ' + str(c['stress'][2,0]) + ' ' + str(c['stress'][2,1]) + ' ' + str(c['stress'][2,2]) + ' ' + ' \n')
        fh.write('c         ' + str(c['c'][0])  + ' ' + str(c['c'][1]) + ' ' + str(c['c'][2]) + ' ' + ' \n')
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
    for c in g.configs['configs']:
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
    efs.make_nl()
  
#@staticmethod
#def efs_add_config():
    
###########################################
#  CLASS e_adjus
###########################################
class e_adjust:

  @staticmethod
  def load():
    if('dft_energy' not in g.inp.keys()):
      return None      
    
    g.dft_energy_adjustments = {}
    
    for k in g.inp['dft_energy'].keys():
    
      atom_label = k
      atom_count = int(g.inp['dft_energy'][k][0])
      relaxed_energy = float(g.inp['dft_energy'][k][1])
      relaxed_energy_unit = str(g.inp['dft_energy'][k][2])
      coh_energy = float(g.inp['dft_energy'][k][3])
      coh_energy_unit = str(g.inp['dft_energy'][k][4])
            
      label_str, label_id = labels.add(atom_label)
           
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
    
    """
# Pd,32,-1264.06398979,Ry,-6.5,eV
    
    g.dft_energy_adjustments = {}
    dft_file = std.path(g.inp['dft']['dir'], g.inp['dft']['e_adjust'])
    if(os.path.isfile(dft_file)): 
# Load csv to list
      fd = std.read_csv(dft_file)
      
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != ''):
          label_str, label_id = labels.add(label_text)
          try:        
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.dft_energy_adjustments[label_id] = {
                                                  'label_id': label_id,
                                                  'label_text': label_str,
                                                  'atom_count': int(row[1]),
                                                  'relaxed_energy': float(row[2]),
                                                  'relaxed_energy_unit': row[3],
                                                  'coh_energy': float(row[4]),
                                                  'coh_energy_unit': row[5],
                                                  'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                  'calc_coh_ev': coh_ev,
                                                  'calc_apaev': apaev,
                                                 }
          except:
            pass  
          
    print()
    print(g.dft_energy_adjustments)
    print()
      
    exit()
    """
      
    """
      g.rd['dft_energy_adjustments'] = {}
      
      dft_energy_adjustments
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != '' and label_text not in g.rd['dft_energy_adjustments'].keys()):
          try:            
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.rd['dft_energy_adjustments'][label_text] = {
                                                      'label': 0,
                                                      'label_text': label_text,
                                                      'atom_count': int(row[1]),
                                                      'relaxed_energy': float(row[2]),
                                                      'relaxed_energy_unit': row[3],
                                                      'coh_energy': float(row[4]),
                                                      'coh_energy_unit': row[5],
                                                      'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                      'calc_coh_ev': coh_ev,
                                                      'calc_apaev': apaev,
                                                     }
          except:
            pass
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
    
    unit_list = [length, energy, force, velocity, pressure, charge_density]
    
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

# There must be an alat value set
        newbp['alat'] = float(bp_inp[k]['alat'])
        newbp['alat'] = units.convert(bp_length, 'ang', newbp['alat'])
                
        try:
          if(bp_inp[k]['type'].lower() == 'sc'):
            newbp['type'] = 1
          elif(bp_inp[k]['type'].lower() == 'bcc'):
            newbp['type'] = 2
          elif(bp_inp[k]['type'].lower() == 'fcc'):
            newbp['type'] = 3
          elif(bp_inp[k]['type'].lower() == 'zb'):
            newbp['type'] = 4
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
          
# Save to list
        g.bulk_properties.append(newbp)

  def make(label_id, label_str):
# b0 bulk modulus
# e0 cohesive energy  (maybe change to ecoh)
# ec elastic constants
#
   
    bp_d ={
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
  
    for bp_n in g.bulk_properties:
#bp_id = bp.add_fcc(6.5, g.bulk_properties[bp_n]['alat'], 1)
#add_bp_config(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, bp_id)
      """
      rcut = g.bulk_properties[bp_n]['rcut']
      alat = g.bulk_properties[bp_n]['alat']
      uv = g.bulk_properties[bp_n]['uv']
      label = g.bulk_properties[bp_n]['label_id']
      type = g.bulk_properties[bp_n]['type']
      expansion = g.bulk_properties[bp_n]['expansion']
      """
      rcut = bp_n['rcut']
      alat = bp_n['alat']
      uv = bp_n['uv']
      label = bp_n['label_id']
      type = bp_n['type']
      expansion = bp_n['expansion']
      
# Add Config
      bp_id = int(bp.add_bp_config(rcut, alat, uv, label, type, expansion))
    
# Add known data
      bp.add_alat(bp_id, bp_n['alat'])
      bp.add_e0(bp_id, bp_n['e0'])
      bp.add_b0(bp_id, bp_n['b0'])
      bp.add_ec(bp_id, bp_n['ec'])
      bp.add_amu_per_crystal(bp_id, bp_n['amu_per_crystal'])
      
      g.bp_ids[bp_id] = {}
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
  def bp_eos_plot():       
  
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

      plt.savefig(g.dirs['eos'] + '/' + 'eos.svg')
      plt.savefig(g.dirs['eos'] + '/' + 'eos.eps')
      
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
               
      plt.savefig(g.dirs['ec'] + '/' + 'ec.svg')
      plt.savefig(g.dirs['ec'] + '/' + 'ec.eps')
      
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
#  CLASS efs_cal
###########################################
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
  
###########################################
#  CLASS bp_cal
###########################################
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
    
    b_props.bp_output()
    b_props.bp_eos_plot()
    
    print('')     
    print('CONFIGS')    
    for n in range(efs.cc):    
      print('Config ' + str(n+1) + ':', efs.config_energy[n,2], efs.energies[n], (efs.config_energy[n,2]-efs.energies[n])**2)
    print('All configs weighted: ' + str(efs.total_rss_weighted))
    
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
    
###########################################
#  CLASS p
###########################################
class pf:
  
  def run():    
    
# Set up EFS and BP modules
    pf.set_up()
    
# Make data structure
    g.pfdata = pf_data.make()
        
# Print start screen
    pf.startup()
    
# Make Pool
    pf_parameters.make_pool()
   
# Run
    pf.run_fit()
    
# Save Results
    pf.save()
    
######################################################
# SET UP - initialise EFS, BP and any other modules
#          and run first calc
######################################################
  def set_up():  
# If a spline fit, convert into a spline with the set number of nodes
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
    bp_calc.get_known()
    
# Rescale density function
    if(g.fit['rescale_density'] == 2):
      rescale_density.run()      
    
  def startup():
    display.clear() 
    
    print("###################################################")      
    print("Starting RSS: ", g.pfdata['rss']['start'])   
    print("###################################################") 
    print("Potential Parameters")
    print("###################################################")   
    potential.print_parameters() 
    
    time.sleep(2.5)

  def run_fit():
  
    for c in range(g.fit['cycles']):
      pf_cycle.run()
  
  def get_rss(top=True):  
    rss = rss_calc.run_calc()  
    
    g.pfdata['rss']['current'] = rss
    g.pfdata['rss']['counter'] += 1
    
    if(rss is not None):    
# Hash of Parameters
      h = pf.param_hash(g.pfdata['params']['current'])    
      
# Store Details
      g.pfdata['rss']['counter_successful'] += 1
      g.pfdata['rss']['since_improvement'] += 1    
      if(g.pfdata['rss']['start'] == None):
        g.pfdata['rss']['start'] = rss
      if(g.pfdata['rss']['best'] == None or rss < g.pfdata['rss']['best']):
        g.pfdata['params']['best'][:] = copy.deepcopy(g.pfdata['params']['current'])
        g.pfdata['best_hash'] = h
        g.pfdata['rss']['best'] = rss
        g.pfdata['bp_best'] = copy.deepcopy(g.bp_results)        
        g.pfdata['rss']['since_improvement'] = 0
        
# Fill in top
      if(top):
        pf.top(g.pfdata['params']['current'], rss)

# DISPLAY
    display.output()    
    return rss
    
  def top(p, rss):  
    g.pfdata['top']['counter'] += 1
    if(g.pfdata['top']['filled']):
      if(rss < g.pfdata['top']['rss'][-1]):
        g.pfdata['top']['rss'][-1] = copy.deepcopy(rss)
        g.pfdata['top']['p'][-1,:] = copy.deepcopy(p)
    else:
      for n in range(g.pfdata['top']['size']):
        breakout = False
        if(g.pfdata['top']['rss'][n] == -1.0):
          g.pfdata['top']['rss'][n] = copy.deepcopy(rss)
          g.pfdata['top']['p'][n,:] = copy.deepcopy(p)
          breakout = True
        if(breakout):
          break
      if(breakout == False):
        g.pfdata['top']['filled'] = True
    
    sort.sort_1d_dp_asc(g.pfdata['top']['rss'][:])
    g.pfdata['top']['rss'] = sort.apply_keytable_1d_dp(g.pfdata['top']['rss'])
    g.pfdata['top']['p'] = sort.apply_keytable_2d_dp(g.pfdata['top']['p'])

#top = {'size': top_size, 'rss': numpy.zeros((top_size,),), 'p': numpy.zeros((top_size, width,),),}
#pass
  
    """
    p = g.pfdata['params']['current']
    for i in range(10):
      if(rss == g.pfdata['top_ten'][i][0]):
        return 0
      if(g.pfdata['top_ten'][i][0] == None):
        g.pfdata['top_ten'][i][0] = rss
        g.pfdata['top_ten'][i][1] = p
        return 0
      if(rss < g.pfdata['top_ten'][i][0]):
        if(i < 9):
          g.pfdata['top_ten'][i+1:9] = copy.deepcopy(g.pfdata['top_ten'][i:8])
        g.pfdata['top_ten'][i][0] = rss
        g.pfdata['top_ten'][i][1] = p
        return 0
    """
    
  def print_top_ten():  
    print()  
    print("=====================")  
    for i in range(10):
      print(i, g.pfdata['top_ten'][i][0])
    print("=====================")  
    print()  
    print()  
    
  def param_hash(p):
    hstr = ''
    for i in range(len(p)):
      hstr = hstr + str(p[i])
    return hashlib.md5(hstr.encode()).hexdigest()
    
  def save():  
   
# Load best
    pf_potential.update(g.pfdata['params']['best'][:])
    
# Display
    g.pfdata['stage'] = 'Finished'
    display.finish()
    
# Output
    potential_output.full()
    
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
    print("# RSS:                ", g.pfdata['rss']['current'], "  [",g.pfdata['rss']['best'],"]") 
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
    display.print_line()
    g.benchmark['configspersec']
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
  
# Central Difference
  def df(p):
    d = numpy.zeros((len(p),),)
    for i in range(len(p)):
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
  def opt(f_rss, p0):
    gd.rss = f_rss   
    gd.last_gamma = 1.0
    p = p0
    
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
    
    return p

###########################################
#  CLASS pf_dat
###########################################
class pf_data:
  
  def make():

    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    width = potential.parameter_count()
    
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
    
    parameters = copy.deepcopy(g.pfdata['params']['start'][:])
    n = 0
    for p in range(pop_size):
      loop = True
      while(loop): 
        n = n + 1
        if(n == 1):
          g.pfdata['params']['pop'][p, :] = g.pfdata['params']['start'][:]
        else:
          g.pfdata['params']['pop'][p, :] = pf_parameters.get_p()
        
# Try - if it fails or rss == None, try next
        try:
# Update
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
# Store
    g.pfdata['params']['current'][:] = p[:]
  
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
      
# Rescale density functions
    if(g.fit['rescale_density'] == 2 and no_rescale == False):
      rescale_density.run()    
    
# Update efs and bp modules
    potential.efs_add_potentials()     # Load potentials
    potential.bp_add_potentials()      # Load potentials

  def take_density(p, p_dens): 
    a = 0
    for fn in range(len(g.pot_functions['functions'])): 
      b = a + g.pot_functions['functions'][fn]['fit_size'] 
      if(g.pot_functions['functions'][fn]['fit_type'] == 1):     # NODE SPLINE 
        pass
      elif(g.pot_functions['functions'][fn]['fit_type'] == 2):   # ANALYTIC          
        if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
          p[a:b] = copy.deepcopy(p_dens[a:b])
      a = b
    return p

###########################################
#  CLASS pf_generatio
###########################################
class pf_generation:

  parents = None
  pop_size = None
  pop_size_half = None
  fresh_size = None

  def run():
  
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
    
#print(g.pfdata['bp'][g.pfdata['best_hash']])
#print(g.bp_known)
      
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
  
    pop_count = g.pfdata['enhance']['top']
    if(pop_count > g.fit['pop_size']):
      pop_count = g.fit['pop_size']
      
#for p in range(pop_count):
#  print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
    
    for p in range(pop_count):
#print(p,gd.rss_out,g.pfdata['params']['pop_rss'][p])
      params = gd.opt(pf_enhance.gd_rss, g.pfdata['params']['pop'][p,:])
      if(gd.rss_out < g.pfdata['params']['pop_rss'][p]):
        g.pfdata['params']['pop_rss'][p] = gd.rss_out
        g.pfdata['params']['pop'][p,:] = numpy.copy(params)
    
        pf_potential.update(g.pfdata['params']['pop'][p,:])
        rss = pf.get_rss(True)
  
#enhance_freq     #enhance_freq
    
  def gd_rss(params):
    pf_potential.update(params)
    rss = pf.get_rss(False)
    return rss 

###########################################
#  CLASS pf_parameter
###########################################
class pf_parameters:

  def get_p():
    pmax = g.pfdata['pool']['pmax']
    pn = g.pfdata['pool']['pn'] % pmax
    p = copy.deepcopy(g.pfdata['pool']['params'][pn, :])                     # Copy parameters
    g.pfdata['pool']['params'][pn, :] = pf_parameters.random_p(p, 0.01)      # Disturb
    g.pfdata['pool']['pn'] = g.pfdata['pool']['pn'] + 1                      # Increment
    return p                                                                 # Return
  
  def make_pool():
    print("Make Pool")
  
    w_start =g.fit['wide_start']
    w_end =g.fit['wide_end']
    pmax = g.fit['pool_size']
    pop_size = g.fit['pop_size']
    fresh_size = g.fit['fresh_size']
    width = potential.parameter_count()
  
    g.pfdata['pool']['params'] = numpy.zeros((pmax,width,),)
    g.pfdata['pool']['pmax'] = pmax
  
# Any Density
    if(g.fit['rescale_density'] == 0):
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(0.0, w)
        g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
    elif(g.fit['rescale_density'] == 1 or g.fit['rescale_density'] == 2):
      sa = g.fit['sane_seeds_a']
      sb = g.fit['sane_seeds_b']
    
      sane_a = numpy.zeros((sa,width,),)
      sane_b = numpy.zeros((sb,width,),)
      
      print("Make Sane A")
      pn = 0
      while(pn < sa):
        params = pf_cycle.random_p()
        pf_potential.update(params, True)
        for fn in range(len(g.pot_functions['functions'])): 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            rho = pf_parameters.estimate_density(fn)
            if( rho > 0.0 and rho <=1.0):
              sane_a[pn,:] = params
              pn = pn + 1
              
      print("Make Sane B")
      pn = 0        
      while(pn < sb):
        r = numpy.random.rand(width)
        params = (1.0 + 0.1 * (0.5-r)) * sane_a[pn%sa,:]
        pf_potential.update(params, True)
        for fn in range(len(g.pot_functions['functions'])): 
          if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
            rho = pf_parameters.estimate_density(fn)
            if( rho > 0.0 and rho <=1.0):
              sane_b[pn,:] = params
              pn = pn + 1
        
      print("Make Pool")
      w = w_start
      w_inc = (w_end - w_start) / (pmax - 1)
      pn = 0
      for n in range(pmax):
        p = pf_parameters.random_p(0.0, w)
        if(g.fit['sane_fraction']>numpy.random.rand()):
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(pf_potential.take_density(p, sane_b[pn%sb,:]))  
        else:
          g.pfdata['pool']['params'][pn, :] = copy.deepcopy(p)
        pn = pn + 1
        w = w + w_inc
        
# Shuffle
    numpy.random.shuffle(g.pfdata['pool']['params'])

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
###########################################
#  MAIN
###########################################
###########################################

def main():

  

# RECORD START TIME
  g.times['start'] = time.time()

  

# SET UP DIRS
  now = datetime.datetime.now()

  date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

  print(date_time)	

  wd_prefix = now.strftime('wd/%Y%m%d_%H%M%S')

  

# Set wd
  g.dirs['wd'] = wd_prefix

  

# Set sub dirs
  if('input' not in g.sub_dirs.keys()):

    g.sub_dirs['input'] = 'input'  

  for sd in g.sub_dirs.keys():

    g.dirs[sd] = g.dirs['wd'] + '/' + g.sub_dirs[sd]

  

# MAKE DIRS
  for d in g.dirs.keys():

    dir = g.dirs[d]

    std.make_dir(dir)

  

# OPEN LOG
  g.log_fh = open(g.dirs['log'] + '/main.log', 'w')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write('\n')

  g.log_fh.write('Script: ' + str(sys.argv[0]) + '\n')

  if(len(sys.argv)>1):

    run_program = False

    try:

      g.inp = read_config.read_file(sys.argv[1])

      g.log_fh.write('Loaded: ' + str(sys.argv[1]) + '\n')

      run_program = True

      

# Copy input
      std.copy(sys.argv[0], g.dirs['input'])

      std.copy(sys.argv[1], g.dirs['input'])

    except:

      g.log_fh.write('Unable to load, exiting\n')

      

# RUN
    if(run_program):

      eampa.run()

    

# CLOSE LOG
  g.times['end'] = time.time()

  g.log_fh.write('\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.write('Duration: ' + str(g.times['end'] - g.times['start']) + '\n')

  g.log_fh.write('###########################################################################\n')

  g.log_fh.close()



# Run
main()

