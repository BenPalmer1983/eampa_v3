class units:
  
  @staticmethod
  def convert(conv_from, conv_to, value_in):
    try:
      value_in = float(value_in)
    except:
      return None
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

    # LENGTH METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

    # ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    }

    # FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
    # VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
    # PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
    # CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }
    
    # TEMPERATURE
    
    
    unit_list = [length, energy, force, velocity, pressure, charge_density]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  
"""  
  @staticmethod
  def convert(conv_from, conv_to, value_in):

    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

    # METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }

    # J
    energy = {
    'J': 1.0,
    'EV': 6.242E18,
    'RY': 4.5874E17,
    }

    if(conv_from in length.keys() and conv_to in length.keys()):
      return round((length[conv_to] / length[conv_from]) * float(value_in),9)

    if(conv_from in energy.keys() and conv_to in energy.keys()):
      return round((energy[conv_to] / energy[conv_from]) * float(value_in),9)
"""

