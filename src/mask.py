######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os
from labels import labels
from pwscf_output import pwscf_output
from potential_functions import potential_functions
from f2py_lib.f_interp import interp

class mask:

  def process():    
    m = None
    for k in g.inp.keys():
      if(k.upper() == "MASK"):
        m = k    
    
    if(m == None):
      return ''
        
    g.mask = {}
    for k in g.inp[m].keys():
      g.mask[k.upper()] = g.inp[m][k].upper()
      
  def get(label):
    label = label.upper()
    if(label in g.mask.keys()):
      return g.mask[label]
    return label
    
    
    
    

