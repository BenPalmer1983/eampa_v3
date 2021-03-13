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
  def cubic_knot_spline_fixed_end(r, p, pf):
    return fnc.cubic_knot_spline_fixed_end_v(r, p, pf)

  @staticmethod
  def make(f, p, pf, xa, xb, node_count):
  
    # Get Nodes
    x = numpy.linspace(xa, xb, node_count)
    y = f(x, p, pf) - 10
    
    print(x)
    print(y)


    xf = numpy.zeros(((len(x) + 4),),)
    xf[0:len(x)] = x[:]
    xf[len(x)] = 6.5
    xf[len(x)+1] = 0.0
    xf[len(x)+2] = 0.0
    xf[len(x)+3] = 1.0

    
    
    plot_x = numpy.linspace(xa, 6.5,101)
    plot_y = f(plot_x, p, pf)
    plot_y_spline = nodes.cubic_knot_spline_fixed_end(plot_x, y, xf)

      


    plt.figure(figsize=(12,8))
    #plt.ylim(-1.0,50.0)
    #plt.plot(plot_x, plot_real_y)
    plt.plot(plot_x, plot_y, 'bx')
    plt.plot(plot_x, plot_y_spline)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('plot1.eps', format='eps')
    plt.close('all') 
   

#TYPE cubic_spline
#PF 2.4 3.2 4.2
#P 11.686859407970 -0.014710740098832 0.47193527075943
p = numpy.asarray([11.686859407970,-0.014710740098832,0.47193527075943])
pf = numpy.asarray([2.4,3.2,4.2])
nodes.make(nodes.cubic_spline, p, pf, 0.0, 4.5, 6)






