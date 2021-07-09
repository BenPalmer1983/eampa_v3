import numpy
import os


class cohesive:

  d = {
      'l_units': None,
      'e_units': None,
      's_units': None,
      'f_units': None,
      'a0_in': None,
      'x_in': [1.0,0.0,0.0],
      'y_in': [0.0,1.0,0.0],
      'z_in': [0.0,0.0,1.0],
      'c_in': [4,4,4],
      'rcut': 6.5,
      'epa': -4.27,
      'coords_in': [],
      'a0': None,
      'x': [1.0,0.0,0.0],
      'y': [0.0,1.0,0.0],
      'z': [0.0,0.0,1.0],
      'c': [1,1,1],
      'coords': [],
      }

  @staticmethod
  def run():
    out_dir = 'out'
    cohesive.make_dir(out_dir)


    fh = open("input.in")
    for line in fh:
      line = cohesive.one_space(line.strip())
      f = line.split(' ')
      if(line != ""):
        if(line[0] == '#'):
          if(line[0:8].upper() == "#L_UNITS"):
            cohesive.d['l_units'] = f[1]
          elif(line[0:8].upper() == "#E_UNITS"):
            cohesive.d['e_units'] = f[1]
          elif(line[0:8].upper() == "#S_UNITS"):
            cohesive.d['s_units'] = f[1]
          elif(line[0:8].upper() == "#F_UNITS"):
            cohesive.d['f_units'] = f[1]
          elif(line[0:5].upper() == "#ALAT"):
            cohesive.d['a0_in'] = float(f[1])
          elif(line[0:2].upper() == "#X"):
            cohesive.d['x_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#Y"):
            cohesive.d['y_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#Z"):
            cohesive.d['z_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#C"):
            cohesive.d['c_in'] = [int(f[1]), int(f[2]), int(f[3])]
        else:
          if(len(f)>= 4):
            cohesive.d['coords_in'].append([f[0],float(f[1]),float(f[2]),float(f[3])])
    fh.close()
    
        # Change coords
    cohesive.d['coords'] = []
    for cx in range(cohesive.d['c_in'][0]):
      for cy in range(cohesive.d['c_in'][1]):
        for cz in range(cohesive.d['c_in'][2]):
          for c in cohesive.d['coords_in']:
            cohesive.d['coords'].append([c[0], (c[1] + cx) / cohesive.d['c_in'][0], (c[2] + cy) / cohesive.d['c_in'][1], (c[3] + cz) / cohesive.d['c_in'][2]])
    

    a0 = float(cohesive.d['c_in'][0]) * float(cohesive.d['a0_in'])
    f = 0.7
    for i in range(401):

      fname = str(i)
      while(len(fname)<4):
        fname = "0" + fname

      fh = open(out_dir + '/' + fname + '.dat', 'w')
      fh.write('#L_UNITS ' + str(cohesive.d['l_units']) + '\n')
      fh.write('#E_UNITS ' + str(cohesive.d['e_units']) + '\n')
      fh.write('#S_UNITS ' + str(cohesive.d['s_units']) + '\n')
      fh.write('#F_UNITS ' + str(cohesive.d['f_units']) + '\n')
      fh.write('#ALAT ' + str(a0 * f) + '\n')
      fh.write('#X ' + str(cohesive.d['x_in'][0]) + ' ' + str(cohesive.d['x_in'][1]) + ' ' + str(cohesive.d['x_in'][2]) + '\n')
      fh.write('#Y ' + str(cohesive.d['y_in'][0]) + ' ' + str(cohesive.d['y_in'][1]) + ' ' + str(cohesive.d['y_in'][2]) + '\n')
      fh.write('#Z ' + str(cohesive.d['z_in'][0]) + ' ' + str(cohesive.d['z_in'][1]) + ' ' + str(cohesive.d['z_in'][2]) + '\n')
      fh.write('#C 1 1 1' + '\n')
      fh.write('#RCUT ' + str(cohesive.d['rcut']) + '\n')
      fh.write('#EPA ' + str(cohesive.d['epa']) + '\n')
      
      for c in cohesive.d['coords']:
        fh.write(c[0] + ' ')
        fh.write(str(c[1]) + ' ')
        fh.write(str(c[2]) + ' ')
        fh.write(str(c[3]) + ' ')
        fh.write('\n')

      fh.close()

      f = f + 0.005


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




cohesive.run()


















