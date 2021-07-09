import numpy
import os


class surface:

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
      'epa': 0.0,
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
    surface.make_dir(out_dir)


    fh = open("input.in")
    for line in fh:
      line = surface.one_space(line.strip())
      f = line.split(' ')
      if(line != ""):
        if(line[0] == '#'):
          if(line[0:8].upper() == "#L_UNITS"):
            surface.d['l_units'] = f[1]
          elif(line[0:8].upper() == "#E_UNITS"):
            surface.d['e_units'] = f[1]
          elif(line[0:8].upper() == "#S_UNITS"):
            surface.d['s_units'] = f[1]
          elif(line[0:8].upper() == "#F_UNITS"):
            surface.d['f_units'] = f[1]
          elif(line[0:5].upper() == "#ALAT"):
            surface.d['a0_in'] = float(f[1])
          elif(line[0:2].upper() == "#X"):
            surface.d['x_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#Y"):
            surface.d['y_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#Z"):
            surface.d['z_in'] = [float(f[1]), float(f[2]), float(f[3])]
          elif(line[0:2].upper() == "#C"):
            surface.d['c_in'] = [int(f[1]), int(f[2]), int(f[3])]
          elif(line[0:4].upper() == "#EPA"):
            surface.d['epa'] = float(f[1])
        else:
          if(len(f)>= 4):
            surface.d['coords_in'].append([f[0],float(f[1]),float(f[2]),float(f[3])])
    fh.close()
    
    # Change coords
    surface.d['coords'] = []
    for cx in range(surface.d['c_in'][0]):
      for cy in range(surface.d['c_in'][1]):
        for cz in range(surface.d['c_in'][2]):
          for c in surface.d['coords_in']:
            surface.d['coords'].append([c[0], (c[1] + cx) / surface.d['c_in'][0], (c[2] + cy) / surface.d['c_in'][1], (c[3] + cz) / surface.d['c_in'][2]])

    # Calculate new UV and a0
    C = numpy.zeros((3,3,),)
    C[0,0] = float(surface.d['c_in'][0])
    C[1,1] = float(surface.d['c_in'][1])
    C[2,2] = float(surface.d['c_in'][2])
    UV = numpy.zeros((3,3,),)
    k = ['x_in', 'y_in', 'z_in']
    for i in range(3):
      for j in range(3):
        UV[i,j] = float(surface.d[k[i]][j])
    M = numpy.matmul(C,UV)
    surface.d['a0'] = surface.d['a0_in'] * M[0,0]
    M = M / M[0,0]


    z_ext = -0.25
    for i in range(12):
      Z = numpy.zeros((3,3,),)
      Z[0,0] = 1.0
      Z[1,1] = 1.0
      Z[2,2] = 1.0 + z_ext
      Z_inv = numpy.linalg.inv(Z)
      UV = numpy.matmul(Z,M)

      fname = str(i)
      while(len(fname)<4):
        fname = "0" + fname

      fh = open(out_dir + '/' + fname + '.dat', 'w')
      fh.write('#L_UNITS ' + str(surface.d['l_units']) + '\n')
      fh.write('#E_UNITS ' + str(surface.d['e_units']) + '\n')
      fh.write('#S_UNITS ' + str(surface.d['s_units']) + '\n')
      fh.write('#F_UNITS ' + str(surface.d['f_units']) + '\n')
      fh.write('#ALAT ' + str(surface.d['a0']) + '\n')
      fh.write('#X ' + str(UV[0,0]) + ' ' + str(UV[0,1]) + ' ' + str(UV[0,2]) + '\n')
      fh.write('#Y ' + str(UV[1,0]) + ' ' + str(UV[1,1]) + ' ' + str(UV[1,2]) + '\n')
      fh.write('#Z ' + str(UV[2,0]) + ' ' + str(UV[2,1]) + ' ' + str(UV[2,2]) + '\n')
      fh.write('#C 1 1 1' + '\n')
      fh.write('#RCUT ' + str(surface.d['rcut']) + '\n')
      fh.write('#EPA ' + str(surface.d['epa']) + '\n')
      
      for c in surface.d['coords']:
        c_coords = numpy.asarray(c[1:])
        c_moved = numpy.matmul(Z_inv, c_coords)


        fh.write(c[0] + ' ')
        fh.write(str(c_moved[0]) + ' ')
        fh.write(str(c_moved[1]) + ' ')
        fh.write(str(c_moved[2]) + ' ')
        fh.write('\n')

      fh.close()

      z_ext = z_ext + 0.05



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




surface.run()


















