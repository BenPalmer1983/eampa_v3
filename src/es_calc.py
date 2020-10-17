######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy
from f2py_lib.f_es import es



class es_calc:



  def run():  
  
    print("Energies: Surface, Vacany etc") 


    # Setup ES
    es.init()
    potential.es_add_potentials()
    es_calc.es_add()
    
    
    #es.energy()
    #for i in range(es.cc):
    #  print(es.config_energy[i,:])
    es.calculate_es()
    

  @staticmethod
  def es_add():  
  
    print(g.es)
    for item in g.es:
      es_id = es.add_es_config(item['rcut'], item['alat'], item['label_id'], item['type'])
  
  
  
  
  
    """ 
    # Test values
    
    rcut = 6.5
    alat = 4.04
    label = 1
    type = 3
      
    # Add Config
    es_id = es.add_es_config(rcut, alat, label, type)
    #print(es_id)
    
    rcut = 6.5
    alat = 3.5
    label = 1
    type = 2
    es_id = es.add_es_config(rcut, alat, label, type)
    #print(es_id)
    """
















  @staticmethod
  def load():
    if('es' not in g.inp.keys()):
      return None
    if('es_file' not in g.inp['es'].keys()):
      return None
    try: 
      dir = g.inp['es']['dir'].strip()
      es_file = std.path(g.inp['es']['dir'], g.inp['es']['es_file'])
    except:
      dir = ""
      es_file = g.inp['es']['es_file']
    
    # Read BP data file
    es_inp = read_config.read_file(es_file)

    
    # Make list
    g.es = []
    
    
    # READ IN UNITS
    try:
      es_pressure = es_inp['units']['pressure']
    except:
      es_pressure = 'GPA'
    try:
      es_length = es_inp['units']['length']
    except:
      es_length = 'ang'
    try:
      es_energy = es_inp['units']['energy']
    except:
      es_energy = 'ev'
      
      
    for k in es_inp.keys():
      if('potlabel' in es_inp[k].keys() and 'alat' in es_inp[k].keys()):
        potlabel = es_inp[k]['potlabel']
        label_str, label_id = labels.add(potlabel)     
        
        new_es = es_calc.make(label_id, label_str)


        # ALAT
        # There must be an alat value set
        new_es['alat'] = units.convert(es_length, 'ang', float(es_inp[k]['alat']))        
          
          
        # UNIT VECTOR  
        new_es['uv'][:,:] = 0.0
        new_es['uv'][0,0] = 1.0
        new_es['uv'][1,1] = 1.0
        new_es['uv'][2,2] = 1.0        
        try: 
          if('uv' in es_inp[k].keys()):
            # CUBIC      
            if(len(es_inp[k]['uv']) == 1):
              new_es['uv'][0,0] = float(es_inp[k]['uv'][0])
              new_es['uv'][1,1] = float(es_inp[k]['uv'][0])
              new_es['uv'][2,2] = float(es_inp[k]['uv'][0])
            # CUBIC      
            if(len(es_inp[k]['uv']) == 3):
              new_es['uv'][0,0] = float(es_inp[k]['uv'][0])
              new_es['uv'][1,1] = float(es_inp[k]['uv'][1])
              new_es['uv'][2,2] = float(es_inp[k]['uv'][2])            
        except:        
          pass
          
          
        # TYPE        
        try:
          if(es_inp[k]['type'].lower() == 'sc'):
            new_es['type'] = 1
          elif(es_inp[k]['type'].lower() == 'bcc'):
            new_es['type'] = 2
          elif(es_inp[k]['type'].lower() == 'fcc'):
            new_es['type'] = 3
          elif(es_inp[k]['type'].lower() == 'zb'):
            new_es['type'] = 4
        except:
          pass 
          
             
        # EXPANSION     
        try: 
          new_es['expansion'] = float(es_inp[k]['expansion'])
        except:
          pass   
        
        
        # RCUT
        try: 
          new_es['rcut'] = units.convert(es_length, 'ang', float(es_inp[k]['rcut']))
        except:
          pass  

        
        # RCUT
        try: 
          new_es['surface_energy'] = units.convert(es_length, 'ang', float(es_inp[k]['surface_energy']))
        except:
          pass 
       
       
       
        # Add to list
        g.es.append(new_es)  

          
          
        
    
    
    
      
      
      
      
      
    #
      
      
      
      
  def make(label_id, label_str):
    # b0 bulk modulus
    # e0 cohesive energy  (maybe change to ecoh)
    # ec elastic constants
    #
   
    es_d ={
           'label_id': label_id,
           'label_str': label_str,
           'alat': None,
           'uv': numpy.zeros((3,3,),),
           'type': 1,
           'expansion': 4,
           'rcut': 6.5,
           'surface_energy': None,
          }
    es_d['uv'][:,:] = 0.0
    es_d['uv'][0,0] = 1.0
    es_d['uv'][1,1] = 1.0
    es_d['uv'][2,2] = 1.0
    return es_d
               
      
      
      