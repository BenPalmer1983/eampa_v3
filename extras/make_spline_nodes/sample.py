from f_fnc import fnc
import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


class nodes:

  @staticmethod
  def cubic_knot_spline(r, p, pf):
    return fnc.cubic_knot_spline_v(r, p, pf)

  @staticmethod
  def cubic_knot_spline_fixed_end(r, p, pf):
    return fnc.cubic_knot_spline_fixed_end_v(r, p, pf)

  @staticmethod
  def make():
    x = numpy.asarray([0.0, 1.0,   2.0,   3.0,  4.0,  5.0,   6.5,0.0,0.0,1.0])
    y = numpy.asarray([1.0, 1.09, -0.015, 0.03, 0.06, 0.039])

    print(nodes.cubic_knot_spline_fixed_end([0.0], y, x))
    print(nodes.cubic_knot_spline_fixed_end([5.0], y, x))
    print(nodes.cubic_knot_spline_fixed_end([6.4], y, x))
    print(nodes.cubic_knot_spline_fixed_end([1.6], y, x))
    print(nodes.cubic_knot_spline_fixed_end([6.49], y, x))
    print(nodes.cubic_knot_spline_fixed_end([6.499], y, x))
    print(nodes.cubic_knot_spline_fixed_end([6.4999999999999], y, x))
    print(nodes.cubic_knot_spline_fixed_end([6.5], y, x))



    plot_x = numpy.linspace(0.0, 6.5,1001)
    plot_y = nodes.cubic_knot_spline_fixed_end(plot_x, y, x)
 
    plt.figure(figsize=(12,8))
    #plt.ylim(-1.0,50.0)
    #plt.plot(plot_x, plot_real_y)
    #plt.plot(plot_x, plot_y)
    plt.plot(x[0:6], y[0:6], 'b+')
    plt.plot(plot_x, plot_y)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('')
    plt.grid(True)
    plt.savefig('plot2.eps', format='eps')
    plt.close('all') 


nodes.make()
