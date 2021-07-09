######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################


class potential:

  fgroup_max = 0
  label_max = 0
            
  pot = {
        'pot_name': '',
        'pot_dir': '',
        'index_file': '',
        'functions': [],
        'p_count': 0,
        'p': None,
        'p_loaded': None,
        'p_var': None,
        'p_lower': None,
        'p_upper': None,
        'oversized': None,
        }
  


  def load():
    print("Load Potential")   


    potential.pot['pot_name'] = ''


 
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

    potential.pot['index_file'] = g.inp['potential']['index_file']
      
    # Read potential index
    potential.read_potential(pot_file)





  # LOAD FROM FILE
  @staticmethod
  def read_potential(file_name):
    potential.pot['pot_dir'] = ''
    if('/' in file_name):
      lst = file_name.split('/')
      for i in range(len(lst) - 1):
        if(i > 0):
          potential.pot['pot_dir'] += '/'
        potential.pot['pot_dir'] += lst[i]
    main.log("Loading: " + str(file_name))        

    index = std.config_file_to_list(file_name)  
    potential.pot['pot_file_name'] = file_name
    pot = potential.pot_function()

    for row in index:    
      if(len(row) > 1 and row[0].upper() == "POTNAME"):
        potential.pot['pot_name'] = row[1]
      elif(len(row) > 1 and row[0].upper() == "F_ON"):
        value = row[1].upper()
        if(value[0].upper() == "N"):
          pot['f_on'] = 0
        if(value[0].upper() == "F"):
          pot['f_on'] = 0
      elif(len(row) > 1 and row[0].upper() == "FILE"):
        pot['file'] = row[1]
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
        else:
          pot['b_text'] = '##ANY##'
          pot['b'] = 0
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
        potential.pot['functions'].append(pot)
        pot = potential.pot_function()


    # Read in files
    for i in range(len(potential.pot['functions'])):
      pf_file = potential.pot['pot_dir'] + '/' + potential.pot['functions'][i]['file']   
      fd = std.config_file_to_list(pf_file)  
      for l in fd:
        if(l[0].upper() == '#TYPE'):
          potential.pot['functions'][i]['f_name'] = l[1].lower().strip()
        elif(l[0].upper().strip() == '#P'):          
          potential.pot['functions'][i]['params'] = numpy.zeros((len(l[1:]),), dtype='float64',)
          for j in range(len(l[1:])):
            potential.pot['functions'][i]['params'][j] = float(l[j+1])
        elif(l[0].upper().strip() == '#PU'):          
          potential.pot['functions'][i]['params_upper'] = numpy.zeros((len(l[1:]),), dtype='float64',)
          for j in range(len(l[1:])):
            potential.pot['functions'][i]['params_upper'][j] = float(l[j+1])
        elif(l[0].upper().strip() == '#PL'):          
          potential.pot['functions'][i]['params_lower'] = numpy.zeros((len(l[1:]),), dtype='float64',)
          for j in range(len(l[1:])):
            potential.pot['functions'][i]['params_lower'][j] = float(l[j+1])
        elif(l[0].upper().strip() == '#PF'):          
          potential.pot['functions'][i]['params_fixed'] = numpy.zeros((len(l[1:]),), dtype='float64',)
          for j in range(len(l[1:])):
            potential.pot['functions'][i]['params_fixed'][j] = float(l[j+1])
        elif(l[0].upper().strip() == '#VR'):       
          try:   
            potential.pot['functions'][i]['vr'] = float(l[1])
          except:
            pass

      if(potential.pot['functions'][i]['vr'] != None and potential.pot['functions'][i]['params_lower'] == None):
        potential.pot['functions'][i]['params_lower'] = potential.pot['functions'][i]['params'] - potential.pot['functions'][i]['vr'] * abs(potential.pot['functions'][i]['params'])
      if(potential.pot['functions'][i]['vr'] != None and potential.pot['functions'][i]['params_upper'] == None):
        potential.pot['functions'][i]['params_upper'] = potential.pot['functions'][i]['params'] + potential.pot['functions'][i]['vr'] * abs(potential.pot['functions'][i]['params'])


    # Integers used by keys
    potential.pot_count = len(potential.pot['functions'])
    potential.label_max = len(g.labels)
    potential.fgroup_max = len(g.fgroups)
    potential.pair_max = potential.key(1, potential.label_max, potential.label_max)
    potential.dens_key_offset = potential.pair_max 
    potential.embe_key_offset = potential.key(2, potential.label_max, potential.label_max, potential.fgroup_max)
    potential.key_max = potential.key(3, potential.label_max, potential.label_max, potential.fgroup_max)

    # Keys
    for i in range(len(potential.pot['functions'])):
      f_type = potential.pot['functions'][i]['f_type_id']
      a = potential.pot['functions'][i]['a']
      b = potential.pot['functions'][i]['b']
      f_group = potential.pot['functions'][i]['f_group']
      potential.pot['functions'][i]['key'] = potential.key(f_type, a, b, f_group)

    # fgroups per label
    potential.f_groups = numpy.zeros((potential.label_max, potential.fgroup_max,), dtype='int32')
    for i in range(len(potential.pot['functions'])):
      if(f_type > 1):
        a = potential.pot['functions'][i]['a'] - 1           # Adjust key Py->Fort
        f_group = potential.pot['functions'][i]['f_group']
        n = 0
        while(True):
          if(potential.f_groups[a,n] == 0):
            potential.f_groups[a,n] = f_group
            break
          elif(potential.f_groups[a,n] == f_group):
            break
          n = n + 1
    #print(potential.f_groups)    


    # Count parameters and make parameter array
    # Only "ON" functions
    potential.pot['p_count'] = 0
    for fn in range(len(potential.pot['functions'])):
      if(potential.pot['functions'][fn]['f_on'] == 1):
        potential.pot['p_count'] = potential.pot['p_count'] + len(potential.pot['functions'][fn]['params'])


    potential.pot['p'] = numpy.zeros((potential.pot['p_count'],), dtype='float64',)
    potential.pot['p_loaded'] = numpy.zeros((potential.pot['p_count'],), dtype='float64',)
    potential.pot['p_var'] = numpy.zeros((potential.pot['p_count'],), dtype='float64',)
    potential.pot['p_lower'] = numpy.zeros((potential.pot['p_count'],), dtype='float64',)
    potential.pot['p_upper'] = numpy.zeros((potential.pot['p_count'],), dtype='float64',)

    a = 0
    for fn in range(len(potential.pot['functions'])):
      if(potential.pot['functions'][fn]['f_on'] == 1):
        b = a + len(potential.pot['functions'][fn]['params'])
        potential.pot['p'][a:b] = numpy.copy(potential.pot['functions'][fn]['params'][:])
        potential.pot['p_loaded'][a:b] = numpy.copy(potential.pot['functions'][fn]['params'][:])
        if(potential.pot['functions'][fn]['params_upper'] is None or potential.pot['functions'][fn]['params_lower'] is None):
          potential.pot['p_var'][a:b] = numpy.zeros((len(potential.pot['p'][a:b]),), dtype='float64',)
          potential.pot['p_lower'][a:b] = numpy.copy(potential.pot['p'][a:b])
          potential.pot['p_upper'][a:b] = numpy.copy(potential.pot['p'][a:b])
        else:
          potential.pot['p_var'][a:b] = numpy.copy((potential.pot['functions'][fn]['params_upper'][:] - potential.pot['functions'][fn]['params_lower'][:]))
          potential.pot['p_lower'][a:b] = numpy.copy(potential.pot['functions'][fn]['params_lower'][:])
          potential.pot['p_upper'][a:b] = numpy.copy(potential.pot['functions'][fn]['params_upper'][:])
        a = b

    potential.pot['oversized'] = None
    if('oversized_parameters' in g.fitting.keys()):
      try:
        potential.pot['oversized'] = numpy.asarray(g.fitting['oversized_parameters'])
      except:
        pass

    potential.print_keys()




  @staticmethod
  def key(f_type, a, b, g=0):
    # PAIR
    if(f_type == 1):
      if(a > potential.label_max or b > potential.label_max):
        return 0
      else:
        min_key = a
        max_key = b
        if(a > b):
          min_key = b
          max_key = a
        return int((max_key*(max_key-1))/2 + min_key)
    # DENS
    elif(f_type == 2):
      if(a > potential.label_max or b > potential.label_max or g > potential.fgroup_max):
        return 0
      else:
        if(b == 0):
          return potential.dens_key_offset + (g-1) * (potential.label_max + potential.pair_max) + a
        elif(b>0):
          min_key = a
          max_key = b
          if(a > b):
            min_key = b
            max_key = a
          return potential.dens_key_offset + (g-1) * (potential.label_max + potential.pair_max) + potential.label_max + int((max_key*(max_key-1))/2 + min_key)
          
          #return potential.dens_key_offset + b + (a - 1) * potential.fgroup_max
    # EMBE
    elif(f_type == 3):
      if(a > potential.label_max or b > potential.label_max or g > potential.fgroup_max):
        return 0
      else:
        return potential.embe_key_offset + (g-1) * (potential.label_max + potential.pair_max) + a
        

  @staticmethod
  def print_keys():
    print("")
    print("DATA")
    print("################################")
    print("Function Count:    ", potential.pot_count)
    print("Label max:         ", potential.label_max)
    print("Fgroup max:        ", potential.fgroup_max)
    print("Pair max:          ", potential.pair_max)
    print("Dens Offset:       ", potential.dens_key_offset)
    print("Embe Offset:       ", potential.embe_key_offset)
    print("Key Max:           ", potential.key_max)

    print("")
    print("Pair")
    print("################################")
    for ka in g.labels.keys():
      kai = g.labels[ka]
      for kb in g.labels.keys():
        kbi = g.labels[kb]
        k = potential.key(1, kai, kbi)
        print(k, ka, kb)
    
    print("")
    print("Density")
    print("################################")
    for kg in g.fgroups.keys():
      kgi = g.fgroups[kg]
      for ka in g.labels.keys():
        kai = g.labels[ka]
        k = potential.key(2, kai, 0, kgi)
        print(k, kg, ka, "Any")
      for ka in g.labels.keys():
        kai = g.labels[ka]
        for kb in g.labels.keys():
          kbi = g.labels[kb]
          k = potential.key(2, kai, kbi, kgi)
          print(k, kg, ka, kb)

    print("")
    print("Embedding")
    print("################################")
    for kg in g.fgroups.keys():
      kgi = g.fgroups[kg]
      for ka in g.labels.keys():
        kai = g.labels[ka]
        k = potential.key(3, kai, 0, kgi)
        print(k, kg, ka)

    print("")
  

  @staticmethod
  def update(p):
    if(potential.pot['p_count'] == len(p)):
      potential.pot['p'] = p
      a = 0
      for fn in range(len(potential.pot['functions'])):
        if(potential.pot['functions'][fn]['f_on'] == 1):
          b = a + len(potential.pot['functions'][fn]['params'])
          potential.pot['functions'][fn]['params'][:] = potential.pot['p'][a:b]
          a = b
      try:
        potential.load_to_efs()
      except:
        pass
      try:
        potential.load_to_bp()
      except:
        pass


  @staticmethod
  def random(p=None, m=None, oversized=False, direction=None):
    # Multiply the range
    if(m == None):
      m = 1.0

    # Randomly pick over sized parameters (defined in input file)
    if(oversized):
      rn = numpy.random.uniform()    
      b = potential.pot['oversized'][0]
      prb = 1.0
      for i in range(1, len(potential.pot['oversized']),1):
        prb = prb * potential.pot['oversized'][i]
        if(rn<prb):
          m = m * b
        else:
          break

    # If none, use the original parameters
    if(p.any() == None):
      p = numpy.copy(potential.pot['p_loaded'])

    # Get random array
    r = 0.5 - numpy.random.rand(potential.pot['p_count'])

    # New parameters
    dp = r * m * potential.pot['p_var']

    # Force direction -1, 0, 1
    try:
      if(direction.any() != None):
        if(len(direction) == len(p)):
          for i in range(len(dp)):
            if(direction[i] == 0.0):
              dp[i] = 0
            elif(direction[i] == 1.0 or direction[i] == -1.0):
              dp[i] = -1.0 * abs(dp[i]) * direction[i]
            elif(direction[i] == 2.0):
              dp[i] = dp[i]
    except:
      pass
    p_new = p + dp

    # Update
    potential.update(p_new)

    # Return
    return p_new
    

  @staticmethod
  def get_p_var():
    return potential.pot['p_var']



  @staticmethod
  def pot_function():
    return {      
    'f_on': 1,  
    'key': 0,
    'file': None,
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=DENS, 3=EMBE
    'f_group': 1,
    'f_group_label': None,
    'f_name': '',
    'params': None,
    'params_fixed': None,
    'vr': None,
    'params_lower': None,
    'params_upper': None,
    }


  






  @staticmethod
  def save(dir_out, p_in=None):

    p = numpy.copy(potential.pot['p'])
    try:
      if(p_in is not None):
        if(len(p_in) == len(p)):
          p = numpy.copy(p_in)
    except:
      pass

    # Make output directory
    std.make_dir(dir_out)
 
    
    fh = open(dir_out +  potential.pot['index_file'], 'w')
    fh.write("POTNAME " + potential.pot['pot_name'] + "\n")
    fh.write("\n")
    for fn in range(len(potential.pot['functions'])):
      if(potential.pot['functions'][fn]['f_on'] == 1):
        fh.write("START\n")
        fh.write("F_ON true\n")
        fh.write("FILE " + potential.pot['functions'][fn]['file'] + "\n")
        if(potential.pot['functions'][fn]['f_type_id'] == 1):
          fh.write("LABEL " + potential.pot['functions'][fn]['a_text'] + " " + potential.pot['functions'][fn]['b_text'] + "\n")
        else:
          fh.write("LABEL " + potential.pot['functions'][fn]['a_text'] + "\n")
        fh.write("F_TYPE " + potential.pot['functions'][fn]['f_type'] + "\n")
        if(potential.pot['functions'][fn]['f_type_id'] > 1):
          fh.write("F_GROUP " + potential.pot['functions'][fn]['f_group_label'] + "\n")
        fh.write("END\n")
        fh.write("\n")
    fh.close()



    a = 0
    for fn in range(len(potential.pot['functions'])):
      if(potential.pot['functions'][fn]['f_on'] == 1):
        b = a + len(potential.pot['functions'][fn]['params'])
        fh = open(dir_out + potential.pot['functions'][fn]['file'], 'w')
        fh.write('#TYPE ' + str(potential.pot['functions'][fn]['f_name']) + '\n') 
        fh.write('#P ')
        for n in range(a,b,1):
          fh.write(str(p[n]) + ' ')
        fh.write('\n')
        fh.write('#PF ')
        for n in range(len(potential.pot['functions'][fn]['params_fixed'])):
          fh.write(str(potential.pot['functions'][fn]['params_fixed'][n]) + ' ')
        fh.write('\n')
        fh.write('#VR ' + str(potential.pot['functions'][fn]['vr']) + '\n') 
        fh.close()
        a = b
    



  @staticmethod
  def plot(dir_out):


    # Make output directory
    std.make_dir(dir_out)
    if(dir_out.strip()[-1] != "/"):
      dir_out = dir_out + "/"

    a = 0
    for fn in range(len(potential.pot['functions'])):
      if(potential.pot['functions'][fn]['f_on'] == 1):
        fname = potential.pot['functions'][fn]['file']
        fname = fname.split('.')
        fname = fname[0]

        b = a + len(potential.pot['functions'][fn]['params'])
        x = numpy.linspace(0.0,10.0,201)
        y = fnc.fv(potential.pot['functions'][fn]['f_name'], x, potential.pot['functions'][fn]['params'], potential.pot['functions'][fn]['params_fixed'])
        a = b

        plt.figure(figsize=(8,6))
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')   
        plt.title(fname)
        plt.plot(x, y, 'k')
        plt.savefig(dir_out + fname + '.eps', type='efs')
        plt.close('all') 
        
        if(potential.pot['functions'][fn]['f_type_id'] == 1):
          plt.figure(figsize=(8,6))
          plt.rc('font', family='serif')
          plt.rc('xtick', labelsize='x-small')
          plt.rc('ytick', labelsize='x-small')   
          plt.plot(x, y, 'k')
          plt.ylim(-5.0,5.0)
          plt.savefig(dir_out + fname + '_zoom.eps', type='efs')
          plt.close('all') 
      



#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

# LOAD INTO MODULES

  @staticmethod
  def load_to_efs():
    #print("Loading")
    efs_potential.pot_count = potential.pot_count
    efs_potential.label_max = potential.label_max
    efs_potential.fgroup_max = potential.fgroup_max
    efs_potential.pair_max = potential.pair_max
    efs_potential.dens_key_offset = potential.dens_key_offset
    efs_potential.embe_key_offset = potential.embe_key_offset
    efs_potential.key_max = potential.key_max
    efs_potential.f_groups[0:potential.label_max,0:potential.fgroup_max] = potential.f_groups[:,:]
    efs_potential.setup()
    #print(efs_potential.f_groups)
    for i in range(potential.pot_count):
      if(potential.pot['functions'][i]['f_on'] == 1):
        key = potential.pot['functions'][i]['key'] - 1   # Adjustment for Python to Fortran
        efs_potential.pot_f_type[key] = potential.pot['functions'][i]['f_type_id']
        efs_potential.pot_option_a[key] = potential.pot['functions'][i]['a']
        efs_potential.pot_option_b[key] = potential.pot['functions'][i]['b']
        efs_potential.pot_f_group[key] = potential.pot['functions'][i]['f_group']
        efs_potential.pot_f_name[key] = potential.pot['functions'][i]['f_name']

        efs_potential.p_count[key] = len(potential.pot['functions'][i]['params'])
        efs_potential.pf_count[key] = len(potential.pot['functions'][i]['params_fixed'])
        pl = len(potential.pot['functions'][i]['params'])
        pfl = len(potential.pot['functions'][i]['params_fixed'])
        efs_potential.params[key,0:pl] = potential.pot['functions'][i]['params'][0:pl]
        efs_potential.paramsfixed[key,0:pfl] = potential.pot['functions'][i]['params_fixed'][0:pfl]
 

        f_type_id = potential.pot['functions'][i]['f_type_id']
        a = potential.pot['functions'][i]['a']
        if(f_type_id == 1):
          b = potential.pot['functions'][i]['b']
        else:
          b = potential.pot['functions'][i]['f_group']


    """
    print("")
    print("Pair")
    print("################################")
    for ka in g.labels.keys():
      kai = g.labels[ka]
      for kb in g.labels.keys():
        kbi = g.labels[kb]
        k = efs_potential.get_pot_key(1, kai, kbi, 0)
        print(k, ka, kb)
    
    print("")
    print("Density")
    print("################################")
    for kg in g.fgroups.keys():
      kgi = g.fgroups[kg]
      for ka in g.labels.keys():
        kai = g.labels[ka]
        k = potential.key(2, kai, 0, kgi)
        print(k, kg, ka, "Any")
      for ka in g.labels.keys():
        kai = g.labels[ka]
        for kb in g.labels.keys():
          kbi = g.labels[kb]
          k = efs_potential.get_pot_key(2, kai, kbi, kgi)
          print(k, kg, ka, kb)

    print("")
    print("Embedding")
    print("################################")
    for kg in g.fgroups.keys():
      kgi = g.fgroups[kg]
      for ka in g.labels.keys():
        kai = g.labels[ka]
        k = potential.key(3, kai, 0, kgi)
        print(k, kg, ka, "Any")
      for ka in g.labels.keys():
        kai = g.labels[ka]
        for kb in g.labels.keys():
          kbi = g.labels[kb]
          k = efs_potential.get_pot_key(3, kai, kbi, kgi)
          print(k, kg, ka, kb)

    print("")
    """

  @staticmethod
  def load_to_bp():
    #print("Loading")
    bp_potential.pot_count = potential.pot_count
    bp_potential.label_max = potential.label_max
    bp_potential.fgroup_max = potential.fgroup_max
    bp_potential.pair_max = potential.pair_max
    bp_potential.dens_key_offset = potential.dens_key_offset
    bp_potential.embe_key_offset = potential.embe_key_offset
    bp_potential.key_max = potential.key_max
    bp_potential.f_groups[0:potential.label_max,0:potential.fgroup_max] = potential.f_groups[:,:]
    bp_potential.setup()
    #print(bp_potential.f_groups)
    for i in range(potential.pot_count):
      if(potential.pot['functions'][i]['f_on'] == 1):
        key = potential.pot['functions'][i]['key'] - 1   # Adjustment for Python to Fortran
        bp_potential.pot_f_type[key] = potential.pot['functions'][i]['f_type_id']
        bp_potential.pot_option_a[key] = potential.pot['functions'][i]['a']
        bp_potential.pot_option_b[key] = potential.pot['functions'][i]['b']
        bp_potential.pot_f_group[key] = potential.pot['functions'][i]['f_group']
        bp_potential.pot_f_name[key] = potential.pot['functions'][i]['f_name']

        bp_potential.p_count[key] = len(potential.pot['functions'][i]['params'])
        bp_potential.pf_count[key] = len(potential.pot['functions'][i]['params_fixed'])
        pl = len(potential.pot['functions'][i]['params'])
        pfl = len(potential.pot['functions'][i]['params_fixed'])
        bp_potential.params[key,0:pl] = potential.pot['functions'][i]['params'][0:pl]
        bp_potential.paramsfixed[key,0:pfl] = potential.pot['functions'][i]['params_fixed'][0:pfl]
 

        f_type_id = potential.pot['functions'][i]['f_type_id']
        a = potential.pot['functions'][i]['a']
        if(f_type_id == 1):
          b = potential.pot['functions'][i]['b']
        else:
          b = potential.pot['functions'][i]['f_group']




