import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter



class embe_fit:



  @staticmethod
  def rss(d, p, f):
    return sum((f(d[:,0], p) - d[:,1])**2)


  @staticmethod
  def embe_a(x, p):
    return p[0] * numpy.sqrt(x) + p[1] * x**2 + p[2] * x**4

  @staticmethod
  def embe_b(x, p):
    return p[0]  + p[1] * x**2 + p[2] * x**4

  @staticmethod
  def embe_c(x, p):
    return p[0] + p[1] * numpy.sqrt(x) + p[2] * x**2 + p[3] * x**4

  @staticmethod
  def read_data(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = embe_fit.one_space(row)
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
  def fit(f, f_size):
    d = embe_fit.read_data('fe_embe.pot')

   
    p = numpy.zeros((f_size,),)
    p_best = p
    rss_best = None
    for n in range(1000000):
      p = p_best + embe_fit.random_p(f_size, 0.9999**n)
      rss = embe_fit.rss(d, p, f)
      if(rss_best == None or rss < rss_best):
        rss_best = rss
        p_best = p
    print(rss_best)

    d_fit = numpy.copy(d)
    d_fit[:,1] = f(d_fit[:,0], p_best)


    plt.figure(figsize=(12,8))
    plt.plot(d[:,0], d[:,1])
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('fit.eps', format='eps')
    plt.close('all') 

    fh = open('fit_embe.pot', "w")
    for i in range(len(d_fit)):
      fh.write(str(d_fit[i,0]) + " " + str(d_fit[i,0]) + "\n")
    fh.close()
 
 
    print(p_best)

embe_fit.fit(embe_fit.embe_c, 4)









