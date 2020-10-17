######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################


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