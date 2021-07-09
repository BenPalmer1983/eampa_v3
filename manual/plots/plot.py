
import numpy
from f_fnc import fnc
import matplotlib.pyplot as plt
import sys
import os



class plot:


  @staticmethod
  def run():

    print("Plot")

    potfiles = []
    for file_name in os.listdir("pots"):
      if(".pot" in file_name):
        potfiles.append(file_name)
     
    for file_name in potfiles:  
      print(file_name)
      
      x = numpy.linspace(0.0, 6.5, 1001)
      y = numpy.zeros((1001,),)
      p_fixed = numpy.asarray([0.0])
      
      fh = open("pots" + "/" + file_name, 'r')
      for line in fh:
        line = plot.one_space(line.strip())
        f = line.split(" ")
        print(f)
        if(f[0] == "#TYPE"):
          f_name = f[1]
        elif(f[0] == "#P"):
          p = numpy.asarray(f[1:])
        elif(f[0] == "#PF"):
          p_fixed = numpy.asarray(f[1:])
        elif(f[0] == "#TITLE"):
          title = line[6:].strip()
        elif(f[0] == "#X"):
          xlabel = line[2:].strip()
        elif(f[0] == "#Y"):
          ylabel = line[2:].strip()  
        elif(f[0] == "#PLOT"):
          plot_name = line[5:].strip()         
      fh.close()
      
      # 
      for n in range(len(x)):
        y[n] = fnc.f(f_name, x[n], p, p_fixed)

      print(title)  


      plt.figure(figsize=(8,6))
      plt.title(title)
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)   
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')
      plt.plot(x, y, color='k', ls='solid')     
      if(max(y) > 1.0e4):
        plt.ylim(-5, 20)  
      if(min(y) < -1.0e4):
        plt.ylim(-50, 50)  
      plt.grid(True)
      plt.savefig("pots/" + plot_name, type='efs')









    exit()

    plot_dir = ""
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

      plot_name = file_name.split(".")
      plot_name = plot_name[0] + '.eps'


      plt.plot(x, y, 'b')
      plt.savefig(plot_dir + "/" + plot_name, type='efs')
      plt.show()

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






