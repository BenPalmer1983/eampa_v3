######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import matplotlib.pyplot as plt
import numpy
import copy
import os
from potential_functions import potential_functions
from potential_output import potential_output
from rescale_density import rescale_density
from f2py_lib.f_interp import interp
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp
from f2py_lib.f_spline import spline


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
      