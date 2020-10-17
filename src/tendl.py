import numpy 
import json
import zlib
import os
import matplotlib.pyplot as plt
from f2py.f_interp import interp

class tendl:

  cache = {}

  @staticmethod
  def read_float(inp):
    out = ''
    if('e' not in inp.lower()):
      for i in range(len(inp)):
        if(i>1 and inp[i] == '+'):
          out = out + 'e'
        elif(i>1 and inp[i] == '-'):
          out = out + 'e-'
        elif(inp[i] != ' '):
          out = out + inp[i]
    else:
      out = inp
    return float(out)


  @staticmethod
  def read_int(inp):
    out = ''
    if('e' not in inp.lower()):
      for i in range(len(inp)):
        if(i>1 and inp[i] == '+'):
          out = out + 'e'
        elif(i>1 and inp[i] == '-'):
          out = out + 'e-'
        elif(inp[i] != ' '):
          out = out + inp[i]
    else:
      out = inp
    return int(numpy.floor(float(out)))

  @staticmethod
  def read_isotope_code(code):
    protons = int(numpy.floor(code/1000))
    nucleons = int(code - 1000 * protons)
    neutrons = nucleons - protons
    return protons, neutrons, nucleons

  @staticmethod
  def mt_change():
    return {
           2: [0,0,0],
           3: [0,0,0],
           4: [0,1,1],
           11: [1,3,4],
           16: [0,2,2],
           17: [0,3,3],
           22: [2,3,5],
           23: [6,7,13],
           24: [2,4,6],
           25: [2,5,7],
           28: [1,1,2],
           29: [2,4,6],
           30: [4,6,10],
           32: [1,2,3],
           33: [1,3,4],
           34: [2,2,4],
           35: [5,6,11],
           36: [5,7,12],
           37: [0,4,4],
           41: [1,2,3],
           42: [1,3,4],
           44: [2,1,3],
           45: [3,3,6],
           }

  @staticmethod
  def read_to_array(block):
    n = 0
    for row in block:
      num = row[0]
      line = row[1] 
      if(num == 3):
        s = tendl.read_int(line[:12])
        #d = numpy.zeros((s,2,),)
        d = []
      elif(num > 3):
        l = []
        l.append(line[0:11])
        l.append(line[11:22])
        l.append(line[22:33])
        l.append(line[33:44])
        l.append(line[44:55])
        l.append(line[55:66])
        m = 0
        while(n<s and m<3):
          #d[n, 0] = read_float(l[2*m])
          #d[n, 1] = read_float(l[2*m+1])
          d.append([tendl.read_float(l[2*m]),tendl.read_float(l[2*m+1])])
          n = n + 1
          m = m + 1
    return d

  @staticmethod
  def make_dir(dir):
    try:
      if(not os.path.exists(dir) and dir.strip() != ''):
        os.mkdir(dir) 
        return True
      return False
    except:
      return False

  @staticmethod
  def convert_file(dir_in, dir_out, file_name, projectile_p, projectile_n):
    f = {}

    fh = open(dir_in + "/" + file_name, 'r')
    line = 0
    mf_last = None
    mt_last = None
    for row in fh:
      mf = int(row[70:72])
      mt = int(row[72:75])
      mat = int(row[66:70])

      if(mt != mt_last or mf != mf_last):
        line = 0
      line = line + 1
      mf_last = mf
      mt_last = mt

      if(mf == 1 and mt == 451 and line == 1):
        this_mat = mat
        code = tendl.read_int(row[0:12])  
        target_protons, target_neutrons, target_nucleons = tendl.read_isotope_code(code)
    
      if(mf not in f.keys()):
        f[mf] = {}

      if(mt not in f[mf].keys()):
        f[mf][mt] = []

      f[mf][mt].append([line, row[:66]])
    fh.close()

    # Read xs
    xs = {}
    mf = 3
    d_mt = tendl.mt_change()

    for mt in d_mt.keys():
      if(mt in f[mf].keys()):
        residual_protons = target_protons - d_mt[mt][0] + projectile_p
        residual_neutrons = target_neutrons - d_mt[mt][1] + projectile_n
        residual_nucleons = residual_protons + residual_neutrons
       
        key = str(target_protons) + "|" + str(target_neutrons) + "|" + str(target_nucleons) + "|" + str(residual_protons) + "|" + str(residual_neutrons) + "|" + str(residual_nucleons)
        if(key in xs.keys()):
          pass
        else:
          dlist = tendl.read_to_array(f[mf][mt])
          xs[key] = {
                  'target_protons': target_protons,
                  'target_neutrons': target_neutrons,
                  'target_nucleons': target_nucleons,
                  'residual_protons': residual_protons,
                  'residual_neutrons': residual_neutrons,
                  'residual_nucleons': residual_nucleons,
                  'data_size': len(dlist),
                  'data': dlist,
                  }
    
    outdir = str(dir_out) + '/' + str(target_protons)
    tendl.make_dir(outdir)
    
    file_name = 'xs_' + str(1000 * projectile_p + (projectile_p + projectile_n)) + '_' + str(1000 * target_protons + target_nucleons) + '.z'
    dat = json.dumps(xs)
    cdat = zlib.compress(dat.encode(), level=9)
    fh = open(outdir + '/' + file_name, 'wb')
    fh.write(cdat)
    fh.close()


  @staticmethod
  def convert_files(dir_in, dir_out, projectile_p, projectile_n):
    tendl.make_dir(dir_out)
    files = os.listdir(dir_in)
    for f in files:
      if('.tendl' in f):
        print(f)
        tendl.convert_file(dir_in, dir_out, f, projectile_p, projectile_n)


  @staticmethod
  def read(dir, pprotons, pneutrons, tprotons, tnucleons):    
    dir = dir + '/' + str(tprotons)
    file_name = 'xs_' + str(1000 * pprotons + (pprotons + pneutrons)) + '_' + str(1000 * tprotons + tnucleons) + '.z'
    fh = open(dir + '/' + file_name, 'rb')
    cdat = ''.encode()
    for r in fh:
      cdat = cdat + r
    d = json.loads(zlib.decompress(cdat).decode())
    for k in d.keys():
      d[k]['data'] = numpy.asarray(d[k]['data'])
    return d

  @staticmethod
  def read_reactions(dir, pprotons, pneutrons, tprotons, tnucleons): 
    r = []
    d = tendl.read(dir, pprotons, pneutrons, tprotons, tnucleons)
    for k in d.keys():
      kl = k.split("|")
      r.append({
                'target_protons': int(kl[0]),
                'target_neutrons': int(kl[1]),
                'target_nucleons': int(kl[2]),
                'residual_protons': int(kl[3]),
                'residual_neutrons': int(kl[4]),
                'residual_nucleons': int(kl[5]),
               })
    return r

  @staticmethod
  def read_reactions_list(dir, pprotons, pneutrons, tprotons, tnucleons): 
    r = []
    d = tendl.read(dir, pprotons, pneutrons, tprotons, tnucleons)
    for k in d.keys():
      kl = k.split("|")
      r.append(kl)
    return r
    
  

  @staticmethod
  def cache_reactions(dir, pprotons, pneutrons, tprotons, tnucleons): 
    rs = tendl.read_reactions_list(dir, pprotons, pneutrons, tprotons, tnucleons)
    for r in rs:
      tendl.read_xs(dir, pprotons, pneutrons, tprotons, tnucleons, r[3], r[5])


  @staticmethod
  def read_xs(dir, pprotons, pneutrons, tprotons, tnucleons, rprotons, rnucleons):
    pprotons = int(pprotons)
    pneutrons = int(pneutrons)
    tprotons = int(tprotons)
    tnucleons = int(tnucleons)
    rprotons = int(rprotons)
    rnucleons = int(rnucleons) 
    code = str(pprotons) + '_' + str(pneutrons) + '_' + str(tprotons) + '_' + str(tnucleons) + '_' + str(rprotons) + '_' + str(rnucleons)
    if(code in tendl.cache.keys()):
      return tendl.cache[code]
    d = tendl.read(dir, pprotons, pneutrons, tprotons, tnucleons) 
    k = str(tprotons) + '|' + str(tnucleons - tprotons) + '|' + str(tnucleons) + '|' + str(rprotons) + '|' + str(rnucleons - rprotons) + '|' + str(rnucleons)
    if(k in d.keys()):      
      tendl.cache[code] = d[k]['data']
      return tendl.cache[code]
      
  @staticmethod
  def plot_xs(dir, pprotons, pneutrons, tprotons, tnucleons, rprotons, rnucleons, plotdir): 
    d = tendl.read_xs(dir, pprotons, pneutrons, tprotons, tnucleons, rprotons, rnucleons)  
    plot_name = str(pprotons) + "_" + str(pneutrons) + "_" + str(tprotons) + "_" + str(tnucleons) + "_" + str(rprotons) + "_" + str(rnucleons)
    plt.clf()
    plt.plot(d[:,0], d[:,1])
    plt.savefig(plotdir + '/' + plot_name + '.eps', format='eps')
    
    
    
  @staticmethod
  def get_xs(dir, pprotons, pneutrons, tprotons, tnucleons, rprotons, rnucleons, energy): 
    d = tendl.read_xs(dir, pprotons, pneutrons, tprotons, tnucleons, rprotons, rnucleons)
    e_min = min(d[:,0])
    e_max = max(d[:,0])    
    if(energy<e_min or energy>e_max):
      return 0.0     
    try:
      xs = float(interp.trap(energy, d[:,0], d[:,1]))
    except: 
      xs = 0.0
    return xs








