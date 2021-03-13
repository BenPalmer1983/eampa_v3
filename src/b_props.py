################################################################
#    Bulk Properties
#
################################################################

import os
import numpy
from std import std
from globals import globals as g
from units import units

"""
Inputs the actual cohesive energy of atoms and relaxed calculated energy so the total energy can be adjusted
"""

class b_props:

  @staticmethod
  def load():
    if('bp' not in g.inp.keys()):
      return None
    if('bp_file' not in g.inp['bp'].keys()):
      return None
    try: 
      dir = g.inp['bp']['dir'].strip()
      bp_file = std.path(g.inp['bp']['dir'], g.inp['bp']['bp_file'])
    except:
      dir = ""
      bp_file = g.inp['bp']['bp_file']

    if(not(os.path.isfile(bp_file))):
      return None
    
    # Read BP data file
    bp_inp = read_config.read_file(bp_file)


    
    # READ IN UNITS
    try:
      bp_pressure = bp_inp['units']['pressure']
    except:
      bp_pressure = 'GPA'
    try:
      bp_length = bp_inp['units']['length']
    except:
      bp_length = 'ang'
    try:
      bp_energy = bp_inp['units']['energy']
    except:
      bp_energy = 'ev'
      
          
    #print(bp_inp)
    
    for k in bp_inp.keys():
      if('potlabel' in bp_inp[k].keys() and 'alat' in bp_inp[k].keys()):
        potlabel = bp_inp[k]['potlabel']
        label_str, label_id = labels.add(potlabel)     
        
        newbp = b_props.make(label_id, label_str)
        newbp['bp_key_str'] = k

        # There must be an alat value set
        newbp['alat'] = float(bp_inp[k]['alat'])
        newbp['alat'] = units.convert(bp_length, 'ang', newbp['alat'])
        #newbp['label_str'] = label_str
                
        try:
          if(bp_inp[k]['type'].lower() == 'sc'):
            newbp['type'] = 1
            newbp['type_text'] = 'Simple Cubic'
          elif(bp_inp[k]['type'].lower() == 'bcc'):
            newbp['type'] = 2
            newbp['type_text'] = 'Body Centered Cubic'
          elif(bp_inp[k]['type'].lower() == 'fcc'):
            newbp['type'] = 3
            newbp['type_text'] = 'Face Centered Cubic'
          elif(bp_inp[k]['type'].lower() == 'zb'):
            newbp['type'] = 4
            newbp['type_text'] = 'Zinc Blende'
        except:
          pass 
          
           
        newbp['uv'][:,:] = 0.0
        newbp['uv'][0,0] = 1.0
        newbp['uv'][1,1] = 1.0
        newbp['uv'][2,2] = 1.0
        
        try: 
          if('uv' in bp_inp[k].keys()):
            # CUBIC      
            if(len(bp_inp[k]['uv']) == 1):
              newbp['uv'][0,0] = float(bp_inp[k]['uv'][0])
              newbp['uv'][1,1] = float(bp_inp[k]['uv'][0])
              newbp['uv'][2,2] = float(bp_inp[k]['uv'][0])
            # CUBIC      
            if(len(bp_inp[k]['uv']) == 3):
              newbp['uv'][0,0] = float(bp_inp[k]['uv'][0])
              newbp['uv'][1,1] = float(bp_inp[k]['uv'][1])
              newbp['uv'][2,2] = float(bp_inp[k]['uv'][2])
            
        except:        
          pass
       
          

          
          
                
        try: 
          newbp['expansion'] = float(bp_inp[k]['expansion'])
        except:
          pass         
        try: 
          newbp['rcut'] = units.convert(bp_length, 'ang', float(bp_inp[k]['rcut']))
        except:
          pass      
        
               
        try:
          newbp['amu_per_crystal'] = float(bp_inp[k]['amu_per_crystal'])
        except:
          pass  
        
        
        try:
          newbp['b0'] = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['b0']))
        except:
          pass  
        try:
          newbp['e0'] = units.convert(bp_energy, 'EV', bp_inp[k]['e0'])
        except:
          pass        
        try:
          if('ec' in bp_inp[k].keys()):
            newbp['ec'] = numpy.zeros((6,6,),)
            
            # CUBIC  C11 C12 C44       
            if(len(bp_inp[k]['ec']) == 3):
              c11 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][0]))
              c12 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][1]))
              c44 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][2]))
                            
              newbp['ec'][0,0] = c11
              newbp['ec'][1,1] = c11
              newbp['ec'][2,2] = c11
              
              newbp['ec'][0,1] = c12
              newbp['ec'][0,2] = c12
              newbp['ec'][1,2] = c12
              newbp['ec'][1,0] = c12
              newbp['ec'][2,0] = c12
              newbp['ec'][2,1] = c12
              
              newbp['ec'][3,3] = c44
              newbp['ec'][4,4] = c44
              newbp['ec'][5,5] = c44
              
            # ORTHORHOMBIC  C11 C22 C33 C44 C55 C66 C12 C13 C23      
            if(len(bp_inp[k]['ec']) == 9):
              c11 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][0]))
              c22 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][1]))
              c33 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][2]))
              c44 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][3]))
              c55 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][4]))
              c66 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][5]))
              c12 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][6]))
              c13 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][7]))
              c23 = units.convert(bp_pressure, 'EV/ANG3', float(bp_inp[k]['ec'][8]))              
              
              newbp['ec'][0,0] = c11
              newbp['ec'][1,1] = c22
              newbp['ec'][2,2] = c33
              
              newbp['ec'][3,3] = c44
              newbp['ec'][4,4] = c55
              newbp['ec'][5,5] = c66
              
              newbp['ec'][0,1] = c12
              newbp['ec'][0,2] = c13
              newbp['ec'][1,2] = c23
              newbp['ec'][1,0] = c12
              newbp['ec'][2,0] = c13
              newbp['ec'][2,1] = c23
              
        except:
          pass
        #print(newbp)  
        # Save to list
        g.bulk_properties.append(newbp)

        
  
  def make(label_id, label_str):
    # b0 bulk modulus
    # e0 cohesive energy  (maybe change to ecoh)
    # ec elastic constants
    #
   
    bp_d ={
           'bp_key_str': '',
           'label_id': label_id,
           'label_str': label_str,
           'alat': None,
           'uv': numpy.zeros((3,3,),),
           'type': 1,
           'expansion': 4,
           'rcut': 6.5,
           'b0': None,
           'e0': None,
           'ec': None,
           'g': None,
           'e': None,
           'poisson': None,
           'amu_per_crystal': None,
           'type': -1,
           'type_text': '',
          }
    bp_d['uv'][:,:] = 0.0
    bp_d['uv'][0,0] = 1.0
    bp_d['uv'][1,1] = 1.0
    bp_d['uv'][2,2] = 1.0
    return bp_d
           
           
      
###########################################################
# F2PY functions
###########################################################
  
  @staticmethod
  def bp_add():  
  
    #print(g.pot_labels)
  
    for bp_i in range(len(g.bulk_properties)):
      bp_n = g.bulk_properties[bp_i]
      label_str = bp_n['label_str'].upper()
      if(label_str in g.pot_labels):

        rcut = bp_n['rcut']
        alat = bp_n['alat']
        uv = bp_n['uv']
        label = bp_n['label_id']
        type = bp_n['type']
        expansion = bp_n['expansion']
      
      
        # Add Config
        bp_id = int(bp.add_bp_config(rcut, alat, uv, label, type, expansion))
        g.bulk_properties[bp_i]['fortran_id'] = bp_id

        # Add known data
        bp.add_alat(bp_id, bp_n['alat'])
        bp.add_e0(bp_id, bp_n['e0'])
        bp.add_b0(bp_id, bp_n['b0'])
        bp.add_ec(bp_id, bp_n['ec'])
        bp.add_amu_per_crystal(bp_id, bp_n['amu_per_crystal'])
      
        g.bp_ids[bp_id] = {}
        g.bp_ids[bp_id]['label_str'] = label_str
        g.bp_ids[bp_id]['rcut'] = rcut
        g.bp_ids[bp_id]['alat'] = alat
        g.bp_ids[bp_id]['uv'] = uv
        g.bp_ids[bp_id]['label'] = label
        g.bp_ids[bp_id]['type'] = type
        g.bp_ids[bp_id]['expansion'] = expansion
        g.bp_ids[bp_id]['e0'] = bp_n['e0']
        g.bp_ids[bp_id]['b0'] = bp_n['b0']
        g.bp_ids[bp_id]['ec'] = bp_n['ec']
        g.bp_ids[bp_id]['amu_per_crystal'] = bp_n['amu_per_crystal']
        g.bp_ids[bp_id]['type'] = bp_n['type']
        g.bp_ids[bp_id]['type_text'] = bp_n['type_text']
        g.bp_ids[bp_id]['bp_key_str'] = bp_n['bp_key_str']
 
    # Add rss multiplication values
    try:
      bp.set_rss_alat(g.inp['rss']['alat'])
    except:
      pass
    try:
      bp.set_rss_e0(g.inp['rss']['e0'])
    except:
      pass
    try:
      bp.set_rss_b0(g.inp['rss']['b0'])
    except:
      pass
    try:
      bp.set_rss_ec(g.inp['rss']['ec'])
    except:
      pass
    try:
      bp.set_rss_g(g.inp['rss']['g'])
    except:
      pass
    try:
      bp.set_rss_e(g.inp['rss']['e'])
    except:
      pass
    try:
      bp.set_rss_v(g.inp['rss']['v'])
    except:
      pass

  
  
  
  @staticmethod
  def bp_output():  
  
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')


      t_pad = 30
      f_pad = 18

      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'LABEL  ' + str(g.bp_ids[bp_id]['label_str'] ), 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Known Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'All units are in ev/Ang unless specified', 1, '', 1)
      std.write_file_line(fh, 'Energy: eV', 1, '', 1)
      std.write_file_line(fh, 'Length: ang', 1, '', 1)
      std.write_file_line(fh, 'Force: eV/ang', 1, '', 1)
      std.write_file_line(fh, 'Pressure: eV/ang3', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.write_file_line(fh, 'Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'Calculated Properties', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Equation of State', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.write_file_line(fh, 'v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.write_file_line(fh, 'e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.write_file_line(fh, 'b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, 'b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Stiffness Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Stiffness (GPA):', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,i,:], f_pad)   
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance Matrix', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Compliance:', t_pad, bp.calc_sc[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Compliance (1/GPA):', t_pad, bp.calc_sc_gpa[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.write_file_line(fh, '', t_pad, bp.calc_sc_gpa[bp_id,i,:], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      
      c11 = 160.230732254e0 * bp.calc_ec[bp_id,0,0]
      c22 = 160.230732254e0 * bp.calc_ec[bp_id,1,1]
      c33 = 160.230732254e0 * bp.calc_ec[bp_id,2,2]
      c44 = 160.230732254e0 * bp.calc_ec[bp_id,3,3]
      c55 = 160.230732254e0 * bp.calc_ec[bp_id,4,4]
      c66 = 160.230732254e0 * bp.calc_ec[bp_id,5,5]
      c12 = 160.230732254e0 * bp.calc_ec[bp_id,0,1]
      c13 = 160.230732254e0 * bp.calc_ec[bp_id,0,2]
      c23 = 160.230732254e0 * bp.calc_ec[bp_id,1,2]
      
      
      
      std.write_file_line(fh, 'Stability', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'C11:', t_pad, c11, f_pad)
      std.write_file_line(fh, 'C11C22 - C12C12:', t_pad, (c11*c22)-(c12*c12), f_pad)
      std.write_file_line(fh, 'C11*C22*C33+2*C12*C13*C23-C11*C23*C23-C33*C12*C12:', t_pad, c11*c22*c33+2*c12*c13*c23-c11*c23*c23-c33*c12*c12, f_pad)
      std.write_file_line(fh, 'C44:', t_pad, c44, f_pad)
      std.write_file_line(fh, 'C55:', t_pad, c55, f_pad)
      std.write_file_line(fh, 'C66:', t_pad, c66, f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
          

      std.write_file_line(fh, 'Bulk Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'b0 reuss:', t_pad, bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight:', t_pad, bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg:', t_pad, bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, 'b0 reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_r[bp_id], f_pad)
      std.write_file_line(fh, 'b0 voight (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_v[bp_id], f_pad)
      std.write_file_line(fh, 'b0 avg (GPA):', t_pad, 160.230732254e0 * bp.calc_b0_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, 'Shear Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'G reuss:', t_pad, bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight:', t_pad, bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg:', t_pad, bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, 'G reuss (GPA):', t_pad, 160.230732254e0 * bp.calc_g_r[bp_id], f_pad)
      std.write_file_line(fh, 'G voight (GPA):', t_pad, 160.230732254e0 * bp.calc_g_v[bp_id], f_pad)
      std.write_file_line(fh, 'G avg (GPA):', t_pad, 160.230732254e0 * bp.calc_g_avg[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1)  
      
      std.write_file_line(fh, 'Young Modulus', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'E:', t_pad, bp.calc_e[bp_id], f_pad)
      std.write_file_line(fh, 'E vec:', t_pad, bp.calc_e_vec[bp_id, :], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Poisson Ratio', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'v:', t_pad, bp.calc_v[bp_id], f_pad)
      std.write_file_line(fh, '', 1, '', 1) 
      
      std.write_file_line(fh, 'Temperatures', 1, '', 1)
      std.write_file_line(fh, '#################', 1, '', 1)
      std.write_file_line(fh, 'Melting (K):', t_pad, bp.calc_melting[bp_id], f_pad)
      #std.write_file_line(fh, 'Debye:', t_pad, bp.calc_debye[bp_id], f_pad)   # Check, calc might be wrong
      std.write_file_line(fh, '', 1, '', 1) 
             
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss], f_pad)
        
       
        
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, 'RSS Weighted', 1, '', 1)
      std.write_file_line(fh, '######################################################', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1) 
      if(bp.known_set[bp_id,0]):
        std.write_file_line(fh, 'Alat rss:', t_pad, bp.rss_w[bp_id,0], f_pad)
      if(bp.known_set[bp_id,1]):
        std.write_file_line(fh, 'e0 rss:', t_pad, bp.rss_w[bp_id,1], f_pad)
      if(bp.known_set[bp_id,2]):
        std.write_file_line(fh, 'b0 rss:', t_pad, bp.rss_w[bp_id,2], f_pad)
      if(bp.known_set[bp_id,3]):
        std.write_file_line(fh, 'Stiffness rss:', t_pad, bp.rss_w[bp_id,3], f_pad)
      if(bp.known_set[bp_id,4]):
        std.write_file_line(fh, 'G rss:', t_pad, bp.rss_w[bp_id,4], f_pad)
      if(bp.known_set[bp_id,5]):
        std.write_file_line(fh, 'E rss:', t_pad, bp.rss_w[bp_id,5], f_pad)
      if(bp.known_set[bp_id,6]):
        std.write_file_line(fh, 'Poisson rss:', t_pad, bp.rss_w[bp_id,6], f_pad)
      std.write_file_line(fh, 'Total rss for ' +str(bp_id+1) + ':', t_pad, bp.rss_total_w[bp_id], f_pad)
      std.write_file_line(fh, 'Total rss for All:', t_pad, [bp.rss_total_rss_w], f_pad)
      
      std.write_file_line(fh, '', 1, '', 1)
      std.write_file_line(fh, '', 1, '', 1)
      
      fh.close()

      
      
  @staticmethod
  def bp_output_terminal():  
  
    t_pad = 30
    f_pad = 18
    
    for bp_id in range(bp.bp_configs_count):
      
      bp_rcut = g.bp_ids[bp_id+1]['rcut']
      bp_alat = g.bp_ids[bp_id+1]['alat']
      bp_uv = g.bp_ids[bp_id+1]['uv']
      bp_label = g.bp_ids[bp_id+1]['label']
      bp_type = g.bp_ids[bp_id+1]['type']
      bp_expansion = g.bp_ids[bp_id+1]['expansion']
      bp_e0 = g.bp_ids[bp_id+1]['e0']
      bp_b0 = g.bp_ids[bp_id+1]['b0']
      bp_ec = g.bp_ids[bp_id+1]['ec']
      bp_amu_per_crystal = g.bp_ids[bp_id+1]['amu_per_crystal']

 


      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)   
      std.print_file_line('############################################################################################################', 1, '', 1)
      std.print_file_line('BP KEY          ' + str(g.bp_ids[bp_id+1]['bp_key_str'].upper()), 1, '', 1)
      std.print_file_line('LABEL           ' + str(g.bp_ids[bp_id+1]['label_str']), 1, '', 1)
      std.print_file_line('TYPE            ' + str(g.bp_ids[bp_id+1]['type_text']), 1, '', 1)
      std.print_file_line('############################################################################################################', 1, '', 1)





      std.print_file_line('', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Known Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('All units are in ev/Ang unless specified', 1, '', 1)
      std.print_file_line('Energy: eV', 1, '', 1)
      std.print_file_line('Length: ang', 1, '', 1)
      std.print_file_line('Force: eV/ang', 1, '', 1)
      std.print_file_line('Pressure: eV/ang3', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
      std.print_file_line('Label (element):', t_pad, bp_label, f_pad)
      std.print_file_line('Structure type:', t_pad, bp_type, f_pad)
      std.print_file_line('Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.print_file_line('Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.print_file_line('AMU per crystal:', t_pad, bp_amu_per_crystal, f_pad)
      std.print_file_line('', 1, '', 1)
      
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
    
      

      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Calculated Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.print_file_line('v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)  
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('', 1, '', 1)
    
    






    
    """ 
    for bp_id in range(bp.bp_configs_count):
      fh = open(g.dirs['results'] + '/' + 'bp_' + str(bp_id + 1) + '.dat', 'w')

      
      
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Known Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('All units are in ev/Ang unless specified', 1, '', 1)
      std.print_file_line('Energy: eV', 1, '', 1)
      std.print_file_line('Length: ang', 1, '', 1)
      std.print_file_line('Force: eV/ang', 1, '', 1)
      std.print_file_line('Pressure: eV/ang3', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Atoms per crystal:', t_pad, bp.known_atoms_per_crystal[bp_id], f_pad)
      std.print_file_line('Expansion:', t_pad, bp.known_expansion[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.known_alat[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.known_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.known_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.known_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.known_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Stiffness (GPA):', t_pad, 160.230732254e0 * bp.known_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, 160.230732254e0 * bp.known_ec[bp_id,i,:], f_pad)
        
      print()
      print()

      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('Calculated Properties', 1, '', 1)
      std.print_file_line('######################################################', 1, '', 1)
      std.print_file_line('', 1, '', 1)
      
      std.print_file_line('Equation of State', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('alat:', t_pad, bp.calc_alat[bp_id], f_pad)
      std.print_file_line('v0:', t_pad, bp.calc_v0[bp_id], f_pad)
      std.print_file_line('e0:', t_pad, bp.calc_e0[bp_id], f_pad)
      std.print_file_line('b0:', t_pad, bp.calc_b0[bp_id], f_pad)
      std.print_file_line('b0/GPA:', t_pad, 160.230732254e0 * bp.calc_b0[bp_id], f_pad)
      std.print_file_line('', 1, '', 1)  
      
      std.print_file_line('Stiffness Matrix', 1, '', 1)
      std.print_file_line('#################', 1, '', 1)
      std.print_file_line('Stiffness:', t_pad, bp.calc_ec[bp_id,0,:], f_pad)
      for i in range(1,6):   
        std.print_file_line('', t_pad, bp.calc_ec[bp_id,i,:], f_pad)
      std.print_file_line('', 1, '', 1)
  """      
        
        

        
  @staticmethod
  def bp_eos_plot(dir):       
  
    
    for bp_id in range(bp.bp_configs_count):
    
      # EQUATION OF STATE
    
      s = bp.calc_sizes[bp_id, 0]
    
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(1, 1, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Equation of State')  
      
      plt.xlabel('Volume (ang3)')    
      plt.ylabel('Energy (eV)')
          
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies[bp_id, 0, 0:s], color='k',  marker="x", ls='')
      plt.plot(bp.calc_volumes[bp_id, 0, 0:s], bp.calc_energies_fit[bp_id, 0, 0:s], color='k', ls='solid')

      
      plt.savefig(dir + '/' + 'eos_' + str(bp_id) + '.svg')
      plt.savefig(dir + '/' + 'eos_' + str(bp_id) + '.eps')
      
      
      # ELASTIC CONSTANTS
      
      plt.clf()
    
      plt.rc('font', family='serif')
      plt.rc('xtick', labelsize='x-small')
      plt.rc('ytick', labelsize='x-small')

      fig, axs = plt.subplots(3, 3, figsize=(12,9))
      fig.tight_layout(pad=5.0)
      fig.suptitle('Elastic Constant Curves')    
    
      for dn in range(9):
        s = bp.calc_sizes[bp_id, dn + 1]
    
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies[bp_id, dn + 1, 0:s],
                                                 color='k',  marker="x", ls='')
        axs[int(numpy.floor(dn/3)), dn % 3].plot(bp.calc_strains[bp_id, dn + 1, 0:s],
                                                 bp.calc_energies_fit[bp_id, dn + 1, 0:s],
                                                 color='k', ls='solid')
        axs[int(numpy.floor(dn/3)), dn % 3].set_title('Distortion D' + str(dn + 1))
        axs[int(numpy.floor(dn/3)), dn % 3].set_xlabel('Strain (Expanded Alat)')
        axs[int(numpy.floor(dn/3)), dn % 3].set_ylabel('Energy (eV)')
               
 
      plt.savefig(dir + '/' + 'ec_' + str(bp_id) + '.svg')
      plt.savefig(dir + '/' + 'ec_' + str(bp_id) + '.eps')
      
    """
    n = globals.d['eos_data_size']
    x = np.linspace(globals.d['eos_data'][1,0], globals.d['eos_data'][1,n-1], 101)    
    y = np.zeros((101,),)
    p = numpy.zeros((4,),)
    p[:] = globals.d['eos_fitting'][:]
    y[:] = run_eos.bm_calc(x[:], p)
    
    
    
    n = globals.d['eos_data_size'] 
    plt.plot(globals.d['eos_data'][1,:n], globals.d['eos_data'][2,:n], color='k',  marker="x", ls='')
    plt.plot(x[:], y[:], color='k', ls='solid')
    """
  
  
  
  