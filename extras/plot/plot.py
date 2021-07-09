
import numpy
from eampa_lib.f_fnc import fnc
import matplotlib.pyplot as plt
import sys



class plot:


  @staticmethod
  def run():

    print("Plot")


    plot_dir = sys.argv[1]
    plot_file = sys.argv[2]


    pot_files = []
    fh = open(plot_dir + "/" + plot_file, 'r')
    for line in fh:
      line = plot.one_space(line.strip())
      f = line.split(" ")
      if(f[0] == "FILE"):
        pot_files.append(f[1])
    fh.close()


    for file_name in pot_files:
      x = numpy.linspace(0.0, 6.5, 1001)
      y = numpy.zeros((1001,),)
      fh = open(plot_dir + "/" + file_name, 'r')
      for line in fh:
        line = plot.one_space(line.strip())
        f = line.split(" ")
        if(f[0] == "#TYPE"):
          f_name = f[1]
        elif(f[0] == "#P"):
          p = numpy.asarray(f[1:])
        elif(f[0] == "#PF"):
          p_fixed = numpy.asarray(f[1:])
      fh.close()

      
      for n in range(len(x)):
        y[n] = fnc.f(f_name, x[n], p, p_fixed)

      plot_name_f = file_name.split(".")
      plot_name = plot_name_f[0] + '.eps'
      plot_name_zoom = plot_name_f[0] + '_zoom.eps'
      
      plt.figure(figsize=(8,6))
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')         
      plt.plot(x, y, 'k')
      plt.savefig(plot_dir + "/" + plot_name, type='efs')
      plt.close()
      
      plt.figure(figsize=(8,6))
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')   
      plt.plot(x, y, 'k')
      plt.ylim(-5.0,5.0)
      plt.savefig(plot_dir + "/" + plot_name_zoom, type='efs')
      plt.close()
      

      

      """
      ymin = min(y)
      ymax = max(y)
      if(ymax > 1.0e8):
        ymax = 10.0
      r = ymax - ymin
      ymin = ymin - 0.05 * ymin
      ymax = ymax + 0.05 * ymax

      
      plt.plot(x, y, 'b')
      plt.savefig(plot_dir + "/" + plot_name, type='efs')
      #plt.show()
      """







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





plot.run()






