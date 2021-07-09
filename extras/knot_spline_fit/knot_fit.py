import sys
import os
import numpy
from eampa_lib.f_fnc import fnc
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import basinhopping



class knot_fit:


  @staticmethod
  def rss(d, p, pf, f):
    return sum((f(d[50:,0], p, pf) - d[50:,1])**2)
    

  @staticmethod
  def rss_opt(p):
    return sum((knot_fit.f(knot_fit.d[:,0], p, knot_fit.pf) - knot_fit.d[:,1])**2)
    
    
  @staticmethod
  def run(file_name, outdir, res=20, d_start=0):
    dir_out = outdir
    knot_fit.make_dir(dir_out)
    d = knot_fit.read_data(file_name)
    
    d = d[d_start:]
    
    min_dx = len(d) // res

    n = []
    fill_gaps = True

    #print(min_dx)
    
    for i in range(len(d)):
      if(i == 0):
        n.append(i)
      elif(i == len(d) - 1):
        n.append(i)
      else:
        ya = d[i-1,1]
        yb = d[i,1]
        yc = d[i+1,1]
        if((ya < yb and yc < yb) or (ya > yb and yc > yb)):
          n.append(i)
      
    if(fill_gaps):
      nn = []
      for i in range(len(n)-1):
        nn.append(n[i])
        if(n[i+1] - n[i] > min_dx): 
          dn = int((n[i+1] - n[i]) / numpy.ceil((n[i+1] - n[i]) / min_dx))
          m = n[i] + dn
          while(m < n[i+1]): 
            nn.append(m)
            m = m + dn
      nn.append(n[-1])
    else:
      nn = n

    p = []
    pf = []
    for i in range(len(nn)):
      pf.append(d[nn[i],0])
      p.append(d[nn[i],1])
      
  
    knot_fit.f = fnc.cubic_knot_spline_v
    knot_fit.d = d[:,:]
    knot_fit.pf = pf    
      

    x = numpy.linspace(d[0,0],d[-1,0],1001)
    y = fnc.cubic_knot_spline_v(x, p, pf)
    
    
    rss = knot_fit.rss_opt(p)
    print(rss)
    
    res = basinhopping(knot_fit.rss_opt, p, niter=1, T=1.0, stepsize=0.5)
    p = res['x']
    rss = knot_fit.rss_opt(p)
    print(rss)
    
    res = minimize(knot_fit.rss_opt, p, method='nelder-mead',
               options={'xatol': 1e-8, 'maxiter': 1000, 'disp': True})
    p = res['x']
    rss = knot_fit.rss_opt(p)
    print(rss)

    res = minimize(knot_fit.rss_opt, p, method='bfgs',
               options={'gtol': 1e-10, 'maxiter': 100, 'disp': True})
    p = res['x']   
    rss = knot_fit.rss_opt(p)
    print(rss)   
     
    

    p_best = p   
    rss_best = rss
              
      


    fh = open(dir_out + "/pot_file.pot", "w")
    fh.write("#TYPE cubic_knot_spline\n")
    fh.write("#P ")
    for i in range(len(p_best)):
      fh.write(str(p_best[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(pf)):
      fh.write(str(pf[i]) + " ")
    fh.write("\n")
    fh.write("#VR 1.0\n")
    fh.close()  
    
    
    
    yb = fnc.cubic_knot_spline_v(x, p_best, pf)
    plt.figure(figsize=(12,8))
    plt.plot(x,y,'r+')
    plt.plot(x,yb,'b+')
    plt.plot(d[:,0],d[:,1])
    plt.grid(True)
    plt.savefig(dir_out + "/plot.eps", format='eps')
    plt.show()
    







  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = knot_fit.one_space(row)
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
  print("python3 knot_fit.py data outdir 20")
  print()
  print("Example:")
  print("python3 knot_fit.py data.pot outdir 20 100")
  exit()
  

filename = str(sys.argv[1]).strip()
outdir = str(sys.argv[2]).strip()
res = int(sys.argv[3])
d_start = 0
if(len(sys.argv) == 5):
  d_start = int(sys.argv[4])
knot_fit.run(filename, outdir, res, d_start)







