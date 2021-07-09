import sys
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from f_fnc import fnc


class test:

  @staticmethod
  def f(x, p, pf):
    return fnc.f(test.fn, x, p, pf)
    
  @staticmethod
  def fv(x, p, pf):
    return fnc.fv(test.fn, x, p, pf)


  @staticmethod
  def run():
   
   

    numpy.random.seed(1)
  
    test.fn = 'cubic_knot_spline_5'
    pf = numpy.asarray([26.0, 26.0, 1.0, 0.0, 6.5, 0.0, 0.0])
    p =  numpy.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 0.2, 0.1, -4.0, -0.2, 0.001, 0,0,0,0,0, 10.0, -5.0]) 

    d_fit = numpy.zeros((1001,2,),)
    d_fit[:,0] = numpy.linspace(0.0, 6.5, 1001)
    d_fit[:,1] = test.fv(d_fit[:,0], p, pf)
    
    plt.figure(figsize=(12,8))
    plt.plot(d_fit[:,0], d_fit[:,1])
    plt.ylim(-5.0,100.0)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    #plt.savefig(plot_file, format='eps')
    plt.show()
    plt.close('all') 
        

test.run()










