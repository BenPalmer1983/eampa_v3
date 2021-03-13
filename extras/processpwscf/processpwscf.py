import os

class processpwscf:


  dft = {}
  configs = {}
  mask = {}


  def run():

    processpwscf.read_input()
    processpwscf.read_files()
    processpwscf.output() 


  def read_input():

    fh = open('dft.in', 'r')
    inp = []
    for line in fh:
      if("=" in line):
        inp.append(line.strip())
    fh.close()  

    for line in inp:
      if("=" in line):
        f = line.split("=")
        label = f[0].upper()
        f = f[1].split(",")
        processpwscf.dft[label] = f

    fh = open('mask.in', 'r')
    inp = []
    for line in fh:
      if("=" in line):
        inp.append(line.strip())
    fh.close()  
   
    for line in inp:
      if("=" in line):
        f = line.split("=")
        label_mask = f[0].upper()
        f = f[1].split(",")
        for fm in f:
          label = fm.upper().strip()
          if(label not in processpwscf.mask):
            processpwscf.mask[label] = label_mask
        
 

  def read_files():
    cwd = os.getcwd() 
    directory = cwd + '/input'
    for filename in os.listdir(directory):
      filepath = os.path.join(directory, filename)
      processpwscf.configs[filename] = processpwscf.read_file(filepath)


  def read_file(filepath):

    config = {
            'nat': None,
            'a0': None,
            'uv': None,
            'coords': None,
            'stress': None,
            'total_energy': None,
            'energy': None,
            }

    f = []
    fh = open(filepath, 'r')
    for line in fh:
      if(line.strip() != ""):
        f.append(line[:-1])

    fh.close()

    n = 0
    while(n<len(f)):
      l = f[n]


      if("lattice parameter" in l):
        config['a0'] = float(l[34:47].strip()) * 0.529 # CONVERT TO ANGSTROM
      elif(config['nat']  == None and "number of atoms/cell" in l):
        config['nat'] = int(l[35:].strip())
        config['coords'] = []
        for n in range(config['nat']):
          config['coords'].append({'label': None, 'x': None, 'y': None, 'z': None, 'fx': None, 'fy': None, 'fz': None})
      elif("crystal axes: (cart. coord. in units of alat)" in l):
        config['uv'] = []
        n = n + 1
        for i in range(3):
          l = f[n]
          config['uv'].append([l[25:36].strip(), l[37:47].strip(), l[47:57].strip()])          
          n = n + 1
      elif("Cartesian axes" in l):
        n = n + 2
        l = f[n]
        i = 0
        while(not "number of k points" in f[n]):
          l = f[n]
          label = l[10:25].strip().upper()
          x = l[40:51].strip()
          y = l[52:63].strip()
          z = l[64:75].strip()
          if(label in processpwscf.mask):
            label = processpwscf.mask[label]
          config['coords'][i]['label'] = label
          config['coords'][i]['x'] = x
          config['coords'][i]['y'] = y
          config['coords'][i]['z'] = z
          i = i + 1
          n = n + 1
      elif("Forces acting on atoms (cartesian axes, Ry/au):" in l):
        n = n + 1
        i = 0
        while("atom" in f[n] and "force" in f[n]):
          fx = float(f[n][34:49].strip()) * 25.710
          fy = float(f[n][50:63].strip()) * 25.710
          fz = float(f[n][64:].strip()) * 25.710
          config['coords'][i]['fx'] = fx
          config['coords'][i]['fy'] = fy
          config['coords'][i]['fz'] = fz
          i = i + 1
          n = n + 1
        
      elif("          total   stress  (Ry/bohr**3) " in l):
        n = n + 1
        config['stress'] = []
        for i in range(3):
          l = f[n]
          sx = float(l[3:15].strip()) * 14583.6
          sy = float(l[15:28].strip()) * 14583.6
          sz = float(l[28:42].strip()) * 14583.6
          config['stress'].append([sx, sy, sz])          
          n = n + 1 
        
      elif("!    total energy   " in l):
        config['total_energy'] = float(l[34:50].strip())  * 13.606  # In EV
 
      n = n + 1


    if(config['total_energy'] is not None):
      coh = 0.0
      relaxed_dft_energy = 0.0
      for i in range(len(config['coords'])):
        e = processpwscf.dft[config['coords'][i]['label']]
        if(e[2].upper() == 'EV'):
          relaxed_dft_energy = relaxed_dft_energy + float(e[1]) / float(e[0])
        elif(e[2].upper() == 'RY'):
          relaxed_dft_energy = relaxed_dft_energy + (float(e[1]) / float(e[0])) * 13.606
        if(e[4].upper() == 'EV'):
          coh = coh + float(e[3])
        elif(e[4].upper() == 'RY'):
          coh = coh + float(e[3]) * 13.606

      config['energy'] = coh + (config['total_energy'] - relaxed_dft_energy)
      config['epa'] = config['energy'] / config['nat']



    return config





  def output():
    print("Output Files")
    
    processpwscf.make_dir('output')
    for f in processpwscf.configs.keys():
      fh = open('output/' + f, 'w')
      print(f)

      a0 = str(processpwscf.configs[f]['a0'])
      uv = processpwscf.configs[f]['uv']
      epa = processpwscf.configs[f]['epa']
      nat = processpwscf.configs[f]['nat']
      coords = processpwscf.configs[f]['coords']
      stress = processpwscf.configs[f]['stress']


      fh.write("/* Converted from PWscf output file */ \n")
      fh.write("#L_UNITS ang \n")
      fh.write("#E_UNITS eV \n")
      fh.write("#S_UNITS GPA \n")
      fh.write("#F_UNITS EV/ANG \n")
      fh.write("#ALAT 	" + processpwscf.float_str(a0)+ "     // a0 \n")
      fh.write("#X 	" + processpwscf.float_str(uv[0][0]) + "  " + processpwscf.float_str(uv[0][1]) + "  " + processpwscf.float_str(uv[0][2]) + "     // X \n")
      fh.write("#Y 	" + processpwscf.float_str(uv[1][0]) + "  " + processpwscf.float_str(uv[1][1]) + "  " + processpwscf.float_str(uv[1][2]) + "     // Y \n")
      fh.write("#Z 	" + processpwscf.float_str(uv[2][0]) + "  " + processpwscf.float_str(uv[2][1]) + "  " + processpwscf.float_str(uv[2][2]) + "     // Z \n")
      fh.write("#SX 	" + processpwscf.float_str(stress[0][0]) + "  " + processpwscf.float_str(stress[0][1]) + "  " + processpwscf.float_str(stress[0][2]) + "     // SX \n")
      fh.write("#SY 	" + processpwscf.float_str(stress[1][0]) + "  " + processpwscf.float_str(stress[1][1]) + "  " + processpwscf.float_str(stress[1][2]) + "     // SY \n")
      fh.write("#SZ 	" + processpwscf.float_str(stress[2][0]) + "  " + processpwscf.float_str(stress[2][1]) + "  " + processpwscf.float_str(stress[2][2]) + "     // SZ \n")


      fh.write("#EPA 	" + processpwscf.float_str(epa) + "     // Energy per Atom \n")
      fh.write("#RCUT 	6.5        // Rcut \n")

      
      for i in range(nat):
        fh.write(coords[i]['label'])
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['x']))
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['y']))
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['z']))
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['fx']))
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['fy']))
        fh.write('   ')
        fh.write(processpwscf.float_str(coords[i]['fz']))
        fh.write('\n')
        

      fh.close()







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
  def float_str(x):
    return str('{:10.5f}'.format(float(x)))



processpwscf.run()