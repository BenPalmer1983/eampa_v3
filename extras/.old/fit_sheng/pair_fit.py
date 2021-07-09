
import numpy
from f_fnc import fnc
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from newtongauss import newtongauss


class pair_fit:


  @staticmethod
  def rss(d, p, f, pf):
    return sum((f(d[:,0], p, pf) - d[:,1])**2)



  @staticmethod
  def ackland_mendelev_pair(r, p, pf):
    return fnc.ackland_mendelev_pair_v(r, p, pf)



  @staticmethod
  def cubic_spline_zbl(r, p, pf):
    return fnc.cubic_spline_zbl_v(r, p, pf)

  @staticmethod
  def cubic_spline_zbl_s(r, p, pf):
    return fnc.cubic_spline_zbl(r, p, pf)
    
    
  @staticmethod
  def pair_spline(r, p, pf):
    return fnc.pair_spline_v(r, p, pf)
    
  @staticmethod
  def pair_spline_s(r, p, pf):
    return fnc.pair_spline(r, p, pf)


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

  @staticmethod
  def random_p(f_size, f=1.0):
    return (f * (0.5-numpy.random.rand(f_size)))

  @staticmethod
  def fit():
    print("Fit")

    f = pair_fit.pair_spline
    fs = pair_fit.pair_spline_s
    d = pair_fit.read_data('fe_pair.pot')
    d_in = numpy.copy(d)
    d = d[320:]

    pf_l = [26,26,1.5,2.0,2.2,2.4,2.7,3.0,3.3,3.7,4.0,4.3,5.0,5.7,6.5]
    pf = numpy.asarray(pf_l)
    f_size = len(pf) - 4
    p = numpy.zeros(f_size)

    
    p_best = p
    rss_best = pair_fit.rss(d, p, f, pf)
    print(rss_best)
    
    for n in range(10000):
      p = pair_fit.random_p(f_size, 0.1)
      rss = pair_fit.rss(d, p, f, pf)
      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p
    
    print(rss_best)


    for n in range(1000):
      p = p_best + pair_fit.random_p(f_size, 0.9999**n)
      rss = pair_fit.rss(d, p, f, pf)
      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p
    print(rss_best)


    
    p_best = newtongauss.fit(f, fs, p_best, pf, d[:,0], d[:,1])
    rss = pair_fit.rss(d, p_best, f, pf)
    print(rss)

    
    d_fit = numpy.copy(d_in)
    d_fit[:,1] = f(d_fit[:,0], p_best, pf)
    print(d_fit)
    
    
   
    plt.figure(figsize=(12,8))
    plt.ylim(-1.0,50.0)
    plt.plot(d_in[:,0], d_in[:,1])
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('fit_pair.eps', format='eps')
    plt.close('all') 
    
    fh = open('fit_pair.pot', "w")
    for i in range(len(d_fit)):
      fh.write(str(d_fit[i,0]) + " " + str(d_fit[i,0]) + "\n")
    fh.close()
 
    fh = open('fit_pair_analytic.pot', "w")
    fh.write("#A\n")
    fh.write("#TYPE pair_spline\n")
    fh.write("#P ")
    for i in range(len(p_best)):
      fh.write(str(p_best[i]) + " ")
    fh.write("\n")
    fh.write("#PF ")
    for i in range(len(pf)):
      fh.write(str(pf[i]) + " ")
    fh.write("\n")
    fh.write("#L 0.0 \n")
    fh.write("#U 7.0 \n")
    fh.close()




#A
#TYPE cubic_spline
#PF 2.4 3.2 4.2
#P 11.686859407970 -0.014710740098832 0.47193527075943
#L 0.0
#U 6.5







pair_fit.fit()

