######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
import os
from f2py_lib.f_fnc import fnc


################################################
#  Define functions here, or define in fnc
#  better for vectorising?
#  Heaviside step etc are in fnc
################################################
 
 
 
class potential_functions:

  ##############################
  #  All the functions have the same format
  #  r which will be a vector/1D array of values
  #  p which will be a vector of parameters
  #  pf which will be a vector of "fixed" parameters - such as r_cut which will be set and fixed, and not varied during the fitting process


  ##############################################
  # PAIR FUNCTIONS
  ##############################################

  # Lennard Jones Potential
  # p[0] = e
  # p[1] = rm
  # f(x) = A * ((B / r)**12 - 2 * (B/r)**6)
  @staticmethod
  def lennard_jones(r, p, pf):
    return fnc.lennard_jones_v(r, p, pf)
    
    
  # Morse Potential
  # p[0] = d
  # p[1] = a
  # p[2] = re
  # f(x) = A * (exp(-2.0D0 * B * (r - C)) - 2.0D0 * exp(-B*(r - C)))
  @staticmethod
  def morse(r, p, pf):
    print(pf)
    return fnc.morse_v(r, p, pf)
    #return 

    
  # Buckingham Potential
  # p[0] = A
  # p[1] = B
  # p[2] = C
  # f(x) = A * exp(-1 * B * r) - C / r**6
  @staticmethod
  def buckingham(r, p, pf):
    return fnc.buckingham_v(r, p, pf)



  ##############################################
  # DENSITY FUNCTIONS
  ##############################################

  # Embedding Finnis-Sinclair
  @staticmethod
  def quadratic_density(r, p, pf):
    return fnc.quadratic_density_v(r, p, pf)
    
    
    
  ##############################################
  # EMBEDDING FUNCTIONS
  ##############################################

  
  # Embedding Finnis-Sinclair
  # f(x) = -A sqrt(rho)
  @staticmethod
  def fs_embedding(r, p, pf):
    return fnc.fs_embedding(r, p, pf)

  
  # Embedding Mendelev
  # f(x) = -sqrt(rho) + A*rho**2
  @staticmethod
  def mendelev_embedding(r, p, pf):
    return fnc.mendelev_embedding_v(r, p, pf)


  # Triple Embedding
  # f(x) = A * sqrt(r) + B * r + C * r**2
  @staticmethod
  def triple_embedding(r, p, pf):
    return fnc.triple_embedding_v(r, p, pf)

  
  # Embedding Ackland (Olsson/Walenius)
  # f(x) = A sqrt(rho) + B rho**2 + C rho**4
  @staticmethod
  def ackland_embedding(r, p, pf):
    return fnc.ackland_embedding_v(r, p, pf)


    



  ##############################################
  # SIMPLE SPLINES
  ##############################################

    
  @staticmethod
  def cubic_spline(r, p, pf):
    return fnc.cubic_spline_v(r, p, pf)
    
  @staticmethod
  def quintic_spline(r, p, pf):
    return fnc.quintic_spline_v(r, p, pf)
    
    
    
    

  ##############################################
  # NODE SPLINES
  ##############################################
   
    
  @staticmethod
  def spline_n_node(r, p, pf):
    return fnc.spline_n_node_v(r, p)
    
    
