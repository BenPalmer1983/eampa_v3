######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os
from potential_functions import potential_functions
from potential import potential
from f2py_lib.f_interp import interp
from f2py_lib.f_efs import efs
from f2py_lib.f_bp import bp

"""
    'f_on': 1,  
    'a_text': '',
    'b_text': '',
    'a': 0,
    'b': 0,
    'f_type': '',             # PAIR, EMBE, DENS
    'f_type_id': 0,           # 1=PAIR, 2=EMBE, 3=DENS
    'f_group': 1,
    'r_cut': 6.5,
    'file': None,
    'function_type': 0,       # 1 tab, 2 analytic
    'f_points': None,         # READ IN TO PYTHON
    'a_type': '',
    'f': None,
    'a_params': None,
    'a_l': 0.0,
    'a_u': 10.0,
    'zoor': 1,
    'points': numpy.zeros((g.tab_size,g.tab_width,),),         # THESE ARE USED BY FORTRAN
"""

class potential_vary:
 
  def vary_all():
    for pn in range(len(g.pot_functions['functions'])):
      if(g.pot_functions['functions'][pn]['function_type'] == 1):
        potential_vary.vary_tabulated(pn)
      elif(g.pot_functions['functions'][pn]['function_type'] == 2):
        potential_vary.vary_analytic(pn)
  
  def vary_tabulated(pn):
    pass
  
  def vary_analytic(pn):
    print(g.pot_functions['functions'][pn])
    
      
      