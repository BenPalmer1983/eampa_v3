
import sys
import numpy
from eampa_lib.f_fnc import fnc
import matplotlib.pyplot as plt



class pair_fit:


  @staticmethod
  def run(file_name, qa, qb, rzbl1, rzbl2, rend, vend, dvdrend, res=20):

    d = pair_fit.read_data(file_name)

    n = 0
    for i in range(len(d)):
      if(d[i, 0] >= rzbl2 and d[i,0] < rend):
        n = n + 1
    
    dd = numpy.zeros((n, 2,),)

    n = 0
    for i in range(len(d)):
      if(d[i, 0] >= rzbl2 and d[i,0] < rend):
        dd[n, 0] = d[i, 0]
        dd[n, 1] = d[i, 1]
        n = n + 1

    n = []

    #print(min_dx)
    
    for i in range(len(dd)):
      if(i == 0):
        n.append(i)
      elif(i == len(dd) - 1):
        n.append(i)
      else:
        ya = dd[i-1,1]
        yb = dd[i,1]
        yc = dd[i+1,1]
        if((ya < yb and yc < yb) or (ya > yb and yc > yb)):
          n.append(i)
      
    if(res>0):
      min_dx = len(dd) // res
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
      pf.append(dd[nn[i],0])
      p.append(dd[nn[i],1])
    pf.append(qa)
    pf.append(qb)
    pf.append(rzbl1)
    pf.append(rend)
    pf.append(vend)
    pf.append(dvdrend)
    
    # qa, qb, rzbl1, rzbl2, rend, vend, dvdrend

    x = numpy.linspace(d[0,0],d[-1,0],1001)
    y = fnc.cubic_knot_spline_fixed_end_pair_v(x, p, pf)

    plt.ylim(-1,10)
    plt.plot(d[:,0],d[:,1])
    plt.plot(dd[:,0],dd[:,1])
    plt.plot(pf[:len(p)],p,'bo')
    plt.plot(x, y,'rx')
    plt.savefig('fit.eps', format='eps')
    plt.close('all') 



    pot_file = file_name.split(".")
    pot_file = pot_file[0] + "_a.pot" 

    fh = open(pot_file, "w")
    fh.write("#TYPE cubic_knot_spline_fixed_end_pair\n")
    fh.write("#P ")
    for i in range(len(p)):
      fh.write(str(p[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(pf)):
      fh.write(str(pf[i]) + " ")
    fh.write("\n")
    fh.write("#VR 1.0\n")
    fh.close()   






  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = pair_fit.one_space(row)
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




if(len(sys.argv) != 4):
  print("How to use:")
  print("python3 fit_pair.py pairfile qa qb")
  print()
  print("Example:")
  print("python3 fit_pair.py al_pair.pot 13 13")
  exit()
  

filename = str(sys.argv[1]).strip()
qa = int(sys.argv[2])
qb = int(sys.argv[2])
pair_fit.run(filename, qa, qb, 0.75, 2.8, 6.5, 0.0, 0.0, 12)






