######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy


# Read pot and config, output these files.

class trial:

  def run():  
  
    print("Trial") 

    print("Saving Configs")
    for cn in range(len(g.configs['configs'])):
      configs.save(cn)
      
    print("Save Potentials")
    potential.save_potential()
      
    potential.plot_python_potentials()
