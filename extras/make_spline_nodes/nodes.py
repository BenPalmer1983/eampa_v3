from f_fnc import fnc
import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


class nodes:

  @staticmethod
  def cubic_spline(r, p, pf):
    return fnc.cubic_spline_v(r, p, pf)
    
  @staticmethod
  def cubic_knot_spline(r, p, pf):
    return fnc.cubic_knot_spline_v(r, p, pf)

  @staticmethod
  def make(f, p, pf, xa, xb, node_count):
    x = numpy.linspace(xa, xb, node_count)
    y = f(x, p, pf)
    
    
    plot_x = numpy.linspace(xa, xb,1001)
    plot_y = f(plot_x, p, pf)


    xf = numpy.zeros(((len(x) + 1),),)
    xf[0:len(x)] = x[:]
    xf[len(x)] = 1.0
    
    plot_x = numpy.linspace(xa, 4.5,100)
    plot_y_spline = nodes.cubic_knot_spline(plot_x, y, xf)

   
    """

    plot_real_y = nodes.cubic_knot_spline(plot_x, y, x)
    """
    plt.figure(figsize=(12,8))
    #plt.ylim(-1.0,50.0)
    #plt.plot(plot_x, plot_real_y)
    plt.plot(plot_x, plot_y_spline)
    plt.plot(x, y)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('plot.eps', format='eps')
    plt.close('all') 
   

#TYPE cubic_spline
#PF 2.4 3.2 4.2
#P 11.686859407970 -0.014710740098832 0.47193527075943
p = numpy.asarray([11.686859407970,-0.014710740098832,0.47193527075943])
pf = numpy.asarray([2.4,3.2,4.2])
nodes.make(nodes.cubic_spline, p, pf, 0.0, 4.5, 11)






