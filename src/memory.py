################################################################
#    Memory Allocation
#
################################################################

import os
from std import std
from globals import globals as g

"""
Inputs the actual cohesive energy of atoms and relaxed calculated energy so the total energy can be adjusted
"""

class memory:

  def run():
  
    
    main.log_title("Memory")
  
    g.memory = {}


    # Defaults
    g.memory['bp'] = {}    
    g.memory['efs'] = {}
      
    
    # For BP module
    try:
      mem = str(g.inp['mem']['bp'])
    except:
      mem = "500MB"
    mem = std.mem_value(mem)    
    g.memory['bp']['mem'] = mem   
    g.memory['bp']['c'] = int(1 * (mem / 7840))
    g.memory['bp']['g'] = int(12 * (mem / 7840))
    g.memory['bp']['nl'] = int(100 * (mem / 7840))
      
    
    # For EFS module
    try:
      mem = str(g.inp['mem']['efs'])
    except:
      mem = "500MB"
    mem = std.mem_value(mem)  
    g.memory['efs']['mem'] = mem    
    g.memory['efs']['c'] = int(1 * (mem / 7840))
    g.memory['efs']['g'] = int(12 * (mem / 7840))
    g.memory['efs']['nl'] = int(100 * (mem / 7840))
    
    
    
    
    main.log("EFS " + str(mem))
    main.log("config size       " + str(g.memory['efs']['c']))
    main.log("ghost size        " + str(g.memory['efs']['g']))
    main.log("nl size           " + str(g.memory['efs']['nl']))

    
    main.log("BP " + str(mem))
    main.log("config size       " + str(g.memory['bp']['c']))
    main.log("ghost size        " + str(g.memory['bp']['g']))
    main.log("nl size           " + str(g.memory['bp']['nl']))





    #print(g.memory)