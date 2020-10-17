########################################################
#  When fitting, rescale density 
#  Max density should be estimated
#  The embedding energy will run from 0.0 to 1.0

from f2py_lib.f_interp import interp
import numpy

# Estimate density function needs to be expanded to include BCC and a variety of lattice parameters

class rescale_density:
  
  def run():
    
    
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        min_rho = min(g.pot_functions['functions'][fn]['points'][:,1])
        g.pot_functions['functions'][fn]['points'][:,1] = g.pot_functions['functions'][fn]['points'][:,1] - min_rho 
        rho = rescale_density.estimate_density(fn)        
        g.pot_functions['functions'][fn]['points'][:,1:3] = (0.5 / rho) * g.pot_functions['functions'][fn]['points'][:,1:3]
        rho = rescale_density.estimate_density(fn)


  def estimate_density(fn):
    r = numpy.zeros((7,),)
    rn = numpy.zeros((7,),)
    r[0] = 7.48332e0
    r[1] = 6.32456e0
    r[2] = 6.92821e0
    r[3] = 4.89898e0
    r[4] = 5.65686e0
    r[5] = 2.82843e0
    r[6] = 4.0e0
    rn[0] = 48
    rn[1] = 24
    rn[2] = 8
    rn[3] = 24
    rn[4] = 12
    rn[5] = 12
    rn[6] = 6    
    rho = 0.0    
    for i in range(7):
      y = interp.search_x(r[i], g.pot_functions['functions'][fn]['points'][:,0], g.pot_functions['functions'][fn]['points'][:,1])
      rho = rho + rn[i] * y
    return rho
    
  def max_densities():    
    density_list = None
    for fn in range(len(g.pot_functions['functions'])): 
      if(g.pot_functions['functions'][fn]['f_type_id'] == 2):
        rho = rescale_density.estimate_density(fn) 
        if(density_list is None):
          density_list = str(rho)
        else:
          density_list = density_list + "/" + str(rho)
    return density_list
        