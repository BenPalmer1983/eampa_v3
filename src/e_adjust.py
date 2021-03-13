################################################################
#    DFT Energy Adjustments
#
################################################################

import os
from std import std
from globals import globals as g
from units import units

"""
Inputs the actual cohesive energy of atoms and relaxed calculated energy so the total energy can be adjusted
"""

class e_adjust:

  @staticmethod
  def load():
    if('dft_energy' not in g.inp.keys()):
      return None      
    
   
    
    g.dft_energy_adjustments = {}
    
    for k in g.inp['dft_energy'].keys():
    
      atom_label = k
      atom_count = int(g.inp['dft_energy'][k][0])
      relaxed_energy = float(g.inp['dft_energy'][k][1])
      relaxed_energy_unit = str(g.inp['dft_energy'][k][2])
      coh_energy = float(g.inp['dft_energy'][k][3])
      coh_energy_unit = str(g.inp['dft_energy'][k][4])
            
      label_str, label_id = labels.add(atom_label)
           
      relaxed_dft_ev = units.convert(relaxed_energy_unit, "EV", relaxed_energy / atom_count) # relaxed per atom in eV
      coh_ev = units.convert(coh_energy_unit, "EV", coh_energy)
      apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
      g.dft_energy_adjustments[label_id] = {
                                            'label_id': label_id,
                                            'label_text': label_str,
                                            'atom_count': atom_count,
                                            'relaxed_energy': relaxed_energy,
                                            'relaxed_energy_unit': relaxed_energy_unit,
                                            'coh_energy': coh_energy,
                                            'coh_energy_unit': coh_energy_unit,
                                            'calc_relaxed_dft_ev': relaxed_dft_ev,
                                            'calc_coh_ev': coh_ev,
                                            'calc_apaev': apaev,
                                           }
    
    print(g.dft_energy_adjustments)

    """
    # Pd,32,-1264.06398979,Ry,-6.5,eV
    
    
    g.dft_energy_adjustments = {}
    dft_file = std.path(g.inp['dft']['dir'], g.inp['dft']['e_adjust'])
    if(os.path.isfile(dft_file)): 
      # Load csv to list
      fd = std.read_csv(dft_file)
      
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != ''):
          label_str, label_id = labels.add(label_text)
          try:        
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.dft_energy_adjustments[label_id] = {
                                                  'label_id': label_id,
                                                  'label_text': label_str,
                                                  'atom_count': int(row[1]),
                                                  'relaxed_energy': float(row[2]),
                                                  'relaxed_energy_unit': row[3],
                                                  'coh_energy': float(row[4]),
                                                  'coh_energy_unit': row[5],
                                                  'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                  'calc_coh_ev': coh_ev,
                                                  'calc_apaev': apaev,
                                                 }
          except:
            pass  
          
    print()
    print(g.dft_energy_adjustments)
    print()
      
    exit()
    """
      
    """
      g.rd['dft_energy_adjustments'] = {}
      
      
      dft_energy_adjustments
      for row in fd:
        label_text = row[0].strip().upper()
        if(label_text != '' and label_text not in g.rd['dft_energy_adjustments'].keys()):
          try:            
            relaxed_dft_ev = units.convert(row[3], "EV", float(row[2]) / int(row[1])) # relaxed per atom in eV
            coh_ev = units.convert(row[5], "EV", float(row[4]))
            apaev = coh_ev - relaxed_dft_ev # Adjustment per atom ev
            g.rd['dft_energy_adjustments'][label_text] = {
                                                      'label': 0,
                                                      'label_text': label_text,
                                                      'atom_count': int(row[1]),
                                                      'relaxed_energy': float(row[2]),
                                                      'relaxed_energy_unit': row[3],
                                                      'coh_energy': float(row[4]),
                                                      'coh_energy_unit': row[5],
                                                      'calc_relaxed_dft_ev': relaxed_dft_ev,
                                                      'calc_coh_ev': coh_ev,
                                                      'calc_apaev': apaev,
                                                     }
          except:
            pass
    """