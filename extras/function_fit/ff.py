import sys
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from eampa_lib.f_fnc import fnc
from scipy.optimize import minimize
from scipy.optimize import basinhopping

class ff:

  @staticmethod
  def rss(d, p, pf, f):
    y = fnc.fv(ff.fn, d[:,0], p, pf)
    return sum((y[:] - d[:,1])**2)


  @staticmethod
  def rss_opt(p):
    y = fnc.fv(ff.fn, ff.d[:,0], p, ff.pf)
    return sum((y[:] - ff.d[:,1])**2)
  

  @staticmethod
  def f(x, p, pf):
    return fnc.f(ff.fn, x, p, pf)
    
  @staticmethod
  def fv(x, p, pf):
    return fnc.fv(ff.fn, x, p, pf)


  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = ff.one_space(row)
        f = row.split(" ")
        
        line = []
        for fn in f:
          line.append(float(fn))
        
        d.append(line)

    d = numpy.asarray(d)
    return d
          

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
  def fit(file_name, funcname, dirout, dstart):
  
    ff.make_dir(dirout)  

    pot_file = file_name.split(".")
    pot_file = dirout + '/' + pot_file[0] + "_a.pot" 
    plot_file = dirout + '/' + pot_file[0] + "_plot.eps" 
    
    ff.d = ff.read_data(file_name)    
    ff.fn = funcname

    if(dstart>0):
      ff.d = ff.d[dstart:]

    if(ff.fn == 'cubic_spline'):
      ff.pf = numpy.asarray([0.5, 1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.5])
      p =  numpy.asarray([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    elif(ff.fn == 'cubic_spline_zbl'):
      ff.pf = numpy.asarray([0.5, 1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.5, 46.0, 46.0,0.2,0.4,1.0])
      p =  numpy.asarray([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    elif(ff.fn == 'cubic_spline_zbl_2'):
      ff.pf = numpy.asarray([1.0, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 6.5, 46.0, 46.0])
      p =  numpy.zeros((len(ff.pf)-2,),)      
    elif(ff.fn == 'cubic_knot_spline_2'):
      ff.pf = numpy.asarray([0.0, 0.5, 0.8, 1.2, 1.5, 2.0, 2.5, 2.8015, 3.029, 3.302, 3.575, 3.848, 4.121, 4.394, 4.667, 4.94, 5.213, 5.486, 5.512, 5.785, 6.058, 6.331, 6.3375, 6.5, 0.0, 0.0, 26.0, 26.0, 1.0])
      p =  numpy.asarray([100.0, 60.0, 15.0, 12.0, 4.0, 2.0, 0.0,  -0.0900170505888, -0.109317871545, -0.10148140479, -0.0903115056083, -0.0817208325078, -0.0742251481766, -0.0566733920941, -0.032521881795, -0.0125538124523, 0.000449171492986, 0.0051164085905, 0.00514853993495, 0.00312301098059, 0.000608339550146, -1.13215997037e-05, -1.1472067638e-05])
      
      ff.pf = numpy.asarray([0.0,0.65,0.923,1.196,1.469,1.495,1.5015,1.547,1.794,2.041,2.288,2.535,2.782,3.029,3.302,3.575,3.848,4.121,4.394,4.667,4.94,5.213,5.486,5.512,5.785,6.058,6.331,6.3375,6.5,0.0,0.0,26.0,26.0,1.0]) 
      p =  numpy.zeros((len(ff.pf)-6,),)
      
    elif(ff.fn == 'slater_4s'):
      ff.pf = numpy.asarray([0.0])
      p =  numpy.asarray([0.0,0.0])
    elif(ff.fn == 'embedding_g'):
      ff.pf = numpy.asarray([0.1])
      p =  numpy.asarray([0.0,0.0])
 
    res = basinhopping(ff.rss_opt, p, niter=10, T=50.0, stepsize=0.5)
    p = res['x']

    res = minimize(ff.rss_opt, p, method='nelder-mead',
               options={'xatol': 1e-8, 'maxiter': 10000, 'disp': True})
    p = res['x']

    res = minimize(ff.rss_opt, p, method='bfgs',
               options={'gtol': 1e-10, 'maxiter': 1000, 'disp': True})
    p = res['x']
    
    d_fit = numpy.copy(ff.d)
    d_fit[:,1] = ff.fv(d_fit[:,0], p, ff.pf)


    plt.figure(figsize=(12,8))
    plt.plot(ff.d[:,0], ff.d[:,1])
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig(plot_file, format='eps')
    plt.show()
    plt.close('all') 
    
    
            
    fh = open(pot_file, "w")
    fh.write("#TYPE " + str(ff.fn) + " \n")
    fh.write("#P ")
    for i in range(len(p)):
      fh.write(str(p[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(ff.pf)):
      fh.write(str(ff.pf[i]) + " ")
    fh.write("\n")
    fh.write("#VR 0.1\n")
    fh.close()  

    
    """
    pf = numpy.asarray([0.5, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.5])
    p =  numpy.asarray([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    

    ff.d = d[100:]
    cubic_spline_fit.pf = pf
    cubic_spline_fit.fn = 'cubic_spline'
    
    
    print(cubic_spline_fit.rss_opt(p))
    res = minimize(cubic_spline_fit.rss_opt, p, method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True})
    p = res['x']
    res = minimize(cubic_spline_fit.rss_opt, p, method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True})
    p = res['x']
    print(cubic_spline_fit.rss_opt(p))
    
    """
    
    
    """
    

    """
    exit()


    
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



    
if(len(sys.argv) < 4):
  print("How to use:")
  print("python3 ff.py datafile function outdir")
  print()
  print()
  print("Example:")
  print("python3 ff.py Pd_plot.den slater4s pddens")
  print("./ff.sh Pd_plot.den slater4s pddens")
  print("./ff.sh Pd_plot.den slater4s pddens 100")
  exit()

filename = str(sys.argv[1]).strip()
funcname = str(sys.argv[2]).strip()
outdir = str(sys.argv[3]).strip()
dstart = 0
if(len(sys.argv)==5):
  dstart = int(sys.argv[4])
ff.fit(filename, funcname, outdir, dstart)









