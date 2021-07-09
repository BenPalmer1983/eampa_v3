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
  
    test.fn = 'cubic_spline'
    pf = numpy.asarray([2.3,2.4,2.5,2.6,2.8,3.0,3.2,3.5,4.0,4.5,5.0,6.5])
    #p =  numpy.asarray([ 102.694439294,-105.249470646,74.0193100636,-45.4016875673, 14.9564151215, 16.1476994746, -20.9977564197, 17.0490125523, -7.07554920275, 1.04087409894, 0.0276703165574 ])
    p =  numpy.asarray([1351.44763875, -772.914944396, -339.549954963, 113.206951842, 2.04448533794, -2.23425865411, -2.86217988307, 1.58917803172, -1.14786166554, 0.468652140934, -0.0130980044832, -0.00480752856994]) 


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










