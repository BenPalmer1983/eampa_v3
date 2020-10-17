################################################################
#    Processing PWscf output file
#
#
#
#
################################################################



import os
import datetime
import re
import sys
import numpy

############################
#  pwscf_output
############################

class pwscf_output:

  def __init__(self, file_in=None):
    self.reset()
    if(file_in != None):
      self.load(file_in)
      self.calcs()

  def reset(self):
  
    self.z = numpy.zeros((3,3))
    
    # Important so store in it's own variable
    self.atom_count = 1   

    # Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "error_code": None,
      "converged": False,
      "converged_in": None,
      
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "species": {},
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
      
      "mass_per_crystal": 0.0,
      "density": [],
    }
    
    # Defaults
    self.xyz_units = 'evang'
    self.stress_units = 'gpa'
    
    # Constants
    self.avogadro = 6.02214086E23

  def scf_settings(self):
    return {
    "bravais_lattice_index": None,
    "alat": None,
    "volume": None,
    "electrons": None, 
    "electrons_up": None, 
    "electrons_down": None, 
    "ecut_wfc": None, 
    "ecut_rho": None, 
    "convergence_threshold": None, 
    "mixing_beta": None, 
    "atomic_species": {},
    }
    
  def scf_crystal(self):
    return {
    "alat": 0.0,
    "cell_parameters": numpy.zeros((3,3)),
    "position_units": None,
    "atomic_labels": [],
    "atomic_positions": numpy.zeros((self.atom_count,3)),    
    "alat_adj": 0.0,
    "cell_parameters_adj": numpy.zeros((3,3)),
    "crystal_positions": numpy.zeros((self.atom_count,3)),
    "cell_volume": 0.0,
    "cell_density": 0.0,
    }

  def scf_results(self):
    return {
    "energy": 0.0,
    "total_force": 0.0,
    "stress": numpy.zeros((3,3)),
    "forces": numpy.zeros((self.atom_count,3)),
    "f_on": False,
    "s_on": False,
    }


  #  Load, and use in another program    
  def load(self, file_name): 
  
    # Load data from file    
    data_file = self.load_from_file(file_name)
    self.d = data_file.split("\n") 
  
    #print(self.d)
  
    # Reset data store
    self.reset()
       
    # Load
    self.load_status()
    self.load_type()
    self.load_times()
    self.load_count()
    self.load_cpuinfo()
    self.load_crystal()
    self.load_scf_settings()
    self.load_results()
    self.load_species()
    
    #print(len(self.data['crystals']))
    #print(len(self.data['results']))
    #self.load_scf('final_scf')
    #self.load_results('initial_scf')
    #self.load_results('final_scf')
    

  def load_status(self):    
    # OK
    ###################################
    self.data['ok'] = False
    counter = 0

    for line in self.d:
      line = line.strip()
      if(pwscf_output.compare(line, "JOB DONE.")):
        self.data['job_done'] = True
      if(pwscf_output.compare(line, "Exit code:")):
        self.data['error'] = True              
        try:
          self.data['error_code'] = int(line[13:-1])
        except:
          pass
      if(pwscf_output.compare(line, "convergence has been achieved")):
        self.data['converged'] = True
        try:
          self.data['converged_in'] = int(line[39:42])
        except:
          pass
    if(self.data['job_done']):
      self.data['ok'] = True
  
  def load_type(self):
    # Calc Type
    ###################################
    self.data['type'] = "SCF"
    for line in self.d:
      line = line.strip()
      if(line[0:23] == "A final scf calculation"):
        self.data['type'] = "VC-RELAX"
    
  def load_times(self):
    # Calc Type
    ###################################
    for line in self.d:
      line = line.strip()
      if(line[0:14] == "PWSCF        :"):
        fa = line.split(':')
        fb = fa[1].split('CPU')
        cpu = fb[0].strip()
        fc = fb[1].split('WALL')
        wall = fc[0].strip()   
        self.data['cpu_time'] = pwscf_output.time_in_seconds(cpu)  
        self.data['wall_time'] = pwscf_output.time_in_seconds(wall)  
        #print(self.data['cpu_time'])
        #print(self.data['wall_time'])
      
  @staticmethod
  def time_in_seconds(inp): 
    h = 0.0
    m = 0.0
    s = 0.0
    if('h' in inp):
      f = inp.split('h')
      try:  
        h = float(f[0])
      except:
        pass
      inp = f[1]
    if('m' in inp):
      f = inp.split('m')
      try:  
        m = float(f[0])
      except:
        pass
      inp = f[1]
    if('s' in inp):
      f = inp.split('s')
      try:  
        s = float(f[0])
      except:
        pass
    return 3600 * h + 60 * m + s
        
    """    
PWSCF        : 15m15.00s CPU    16m33.12s WALL 
PWSCF        :    28.10s CPU        29.31s WALL 
PWSCF        :     1h29m CPU        1h32m WALL         
    """
   
  
  def load_count(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "number of atoms/cell      =")):
        count = pwscf_output.extract(line, "=", "", "i")  
        try:
          self.atom_count = int(count)
          return self.atom_count 
        except:
          return 0
          
  
  def load_cpuinfo(self):
    counter = 0
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
          elif(pwscf_output.compare(line, "Threads/MPI process:")):
            self.data['threads_per_mpi_process'] = pwscf_output.extract(line, ":", "", "i") 
   
   
  def load_species(self):
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(line[0:14] == "atomic species"):
        n, line, line_uc = self.next_line(n, self.d)
        while(line.strip() != ""):      
          line = pwscf_output.single_spaces(line).strip()
          f = line.split(" ")
          self.data['scf_settings'][f[0]] = [float(f[1]), float(f[2])]   # Valence electrons, atomic mass
          n, line, line_uc = self.next_line(n, self.d)
        break
      
      
      #self.data['scf_settings']
      #
      #  atomic species   valence    mass     pseudopotential
      #  Al             3.00    26.98200     Al( 1.00)
 
   
  def load_scf_settings(self):
    self.data['scf_settings'] = self.scf_settings()  
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      if(pwscf_output.compare(line, "bravais-lattice index     =")):
        self.data['scf_settings']['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
      elif(pwscf_output.compare(line, "lattice parameter (alat)  =")):
        self.data['scf_settings']['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
      elif(pwscf_output.compare(line, "unit-cell volume          =")):
        self.data['scf_settings']['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
      elif(pwscf_output.compare(line, "number of atoms/cell      =")):
        self.data['scf_settings']['nat'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of atomic types    =")):
        self.data['scf_settings']['types'] = pwscf_output.extract(line, "=", "", "i")  
      elif(pwscf_output.compare(line, "number of electrons       =")):
        str_e = pwscf_output.extract(line, "=", "", "s")
        e, eu, ed = pwscf_output.electron_string(str_e)
        self.data['scf_settings']['electrons'] = e
        self.data['scf_settings']['electrons_up'] = eu
        self.data['scf_settings']['electrons_down'] = ed
      elif(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
        self.data['scf_settings']['ecut_wfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "charge density cutoff     =")):
        self.data['scf_settings']['ecut_rho'] = pwscf_output.extract(line, "=", "Ry", "f")  
      elif(pwscf_output.compare(line, "convergence threshold     =")):
        self.data['scf_settings']['convergence_threshold'] = pwscf_output.extract(line, "=", "", "f")  
      elif(pwscf_output.compare(line, "mixing beta               =")):
        self.data['scf_settings']['mixing_beta'] = pwscf_output.extract(line, "=", "", "f")  
      elif(("atomic species" in line) and ("valence" in line) and ("mass" in line) and ("pseudopotential" in line)):
        loop = True
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)
          if(line.strip() == ""):
            loop = False
          else:  
            line = pwscf_output.single_spaces(line).strip()
            line_arr = line.split(" ")
            self.data['scf_settings']['atomic_species'][line_arr[0]] = {}
            self.data['scf_settings']['atomic_species'][line_arr[0]]['valence'] = line_arr[1]
            self.data['scf_settings']['atomic_species'][line_arr[0]]['mass'] = line_arr[2]


        # End of file/loop
        n = len(self.d)


  ###################################
  # LOAD CRYSTALS FROM OUTPUT FILE
  ###################################

  def load_crystal(self):
    # Make new list for crystals
    self.data['crystals'] = []
    
    # FIRST
    n = 0
    crystal = self.scf_crystal()
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)  
      if(line[0:10] == "celldm(1)="):
        crystal['alat'] = float(line[10:21].strip())
      elif(pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(pwscf_output.compare(line.strip(), "Cartesian axes")):
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
        # Unit
        crystal['position_units'] = "alat"
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
            # Increment
            k = k + 1
        
        # Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
        
        # Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
    
    # MIDDLE
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:15] == "CELL_PARAMETERS"):
        # Create
        crystal = self.scf_crystal()
        
        # Get alat
        line_arr = line.split("=")
        line_arr = line_arr[1].split(")")
        crystal['alat'] = float(line_arr[0].strip())
        
        #Cell Parameters
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          line = pwscf_output.single_spaces(line)
          fields = line.split(" ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
      elif(line[0:9] == "density ="):  
        crystal['density'] = float(line[10:23])      
            
      elif(line[0:16] == "ATOMIC_POSITIONS"):     
        
        # Unit
        crystal['position_units'] = "crystal"
        
        # Read Coords
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          elif(line.strip() == "End final coordinates"):
            loop = False
          else:
            line = pwscf_output.single_spaces(line)
            line_arr = line.split(" ")
            crystal['atomic_labels'].append(line_arr[0])
            
            crystal['atomic_positions'][k, 0] = float(line_arr[1])
            crystal['atomic_positions'][k, 1] = float(line_arr[2])
            crystal['atomic_positions'][k, 2] = float(line_arr[3]) 
            
            # Increment
            k = k + 1
            
        # Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
        # Add/Save
        self.data['crystals'].append(crystal)
    
    
    # END
    n = 0
    crystal = self.scf_crystal()
    d = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)   
      if(line[0:10] == "celldm(1)="):
        d = d + 1
        if(d == 2):
          crystal['alat'] = float(line[10:21].strip())
      elif(d == 2 and pwscf_output.compare(line.strip(), "crystal axes:")):
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)
          fields = pwscf_output.extract(line, "= (", ")", "s", " ")
          for i in range(len(fields)):
            crystal['cell_parameters'][j, i] = float(fields[i])
     
      elif(d == 2 and pwscf_output.compare(line.strip(), "Cartesian axes")):      
        
        # Unit
        crystal['position_units'] = "alat"
        
        # Read coords
        n, line, line_uc = self.next_line(n, self.d)
        n, line, line_uc = self.next_line(n, self.d)
        
        loop = True 
        k = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)   
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("tau(")
            label = line_arr[0][-15:].strip()
            crystal['atomic_labels'].append(label)
            
            coords = line_arr[1]
            x = float(coords[9:21])
            y = float(coords[22:33])
            z = float(coords[34:44])
            crystal['atomic_positions'][k, 0] = x
            crystal['atomic_positions'][k, 1] = y
            crystal['atomic_positions'][k, 2] = z 
            
            # Increment
            k = k + 1
            
            
        # Cell Volume
        crystal['cell_volume'] = pwscf_output.cell_volume(crystal['alat'], crystal['cell_parameters'][:,:])    
            
        # Add/Save
        self.data['crystals'].append(crystal)
        n = len(self.d)
    
        
    # Loop through crystals
    for i in range(len(self.data['crystals'])):
      # Adjust alat and cell_parameters so celldm(1)=1.0
      factor = 1.0 / self.data['crystals'][i]['cell_parameters'][0, 0]
      
      self.data['crystals'][i]['alat_adj'] = self.data['crystals'][i]['alat'] * self.data['crystals'][i]['cell_parameters'][0, 0]
      self.data['crystals'][i]['cell_parameters_adj'][:, :] = factor * self.data['crystals'][i]['cell_parameters'][:, :]
    
    
      # Make crystal_positions
      if(self.data['crystals'][i]['position_units'] == 'crystal'):
        self.data['crystals'][i]['crystal_positions'][:,:] = self.data['crystals'][i]['atomic_positions'][:,:] 
      elif(self.data['crystals'][i]['position_units'] == 'alat'):  
        minv = numpy.linalg.inv(self.data['crystals'][i]['cell_parameters'][:, :])
        for j in range(len(self.data['crystals'][i]['atomic_positions'])):
          self.data['crystals'][i]['crystal_positions'][j, :] = numpy.matmul(minv[:,:], self.data['crystals'][i]['atomic_positions'][j, :])



  def load_results(self):
  
    # Make new list for results
    self.data['results'] = []
    
    n = 0
    while(n < len(self.d)):
      n, line, line_uc = self.next_line(n, self.d)
      
      # READ ENERGY
      if(pwscf_output.compare(line, "!    total energy")):
        # Create dictionary
        results = self.scf_results()  
        results['energy'] = pwscf_output.extract(line, "=", "Ry", "f") 
        
      # READ FORCES  
      elif(pwscf_output.compare(line, "Forces acting on atoms")):
        n, line, line_uc = self.next_line(n, self.d)  
        loop = True
        f = 0
        while(loop):
          n, line, line_uc = self.next_line(n, self.d)  
          if(line.strip() == ""):
            loop = False
          else:
            line_arr = line.split("force =")
            fields = pwscf_output.single_spaces(line_arr[1].strip()).split(" ")
            results['forces'][f,0] = float(fields[0])
            results['forces'][f,1] = float(fields[1])
            results['forces'][f,2] = float(fields[2])
            f = f + 1
        if(f>0):
          results['f_on'] = True
        
      # READ TOTAL FORCE 
      elif(pwscf_output.compare(line, "Total force =")):
        results['total_force'] = pwscf_output.extract(line, "=", "T", "f")
        
      # READ STRESS
      elif(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):  
        results['s_on'] = True      
        for j in range(3):              
          n, line, line_uc = self.next_line(n, self.d)  
          fields = pwscf_output.extract(line, "", "", "f", " ", True)  
          results['stress'][j,0] = fields[3] 
          results['stress'][j,1] = fields[4] 
          results['stress'][j,2] = fields[5]

        #SAVE
        self.data['results'].append(results)


             
  def calcs(self):
  
    # MASS PER CELL
    self.data['mass_per_crystal'] = 0.0    
    for l in self.data['crystals'][0]['atomic_labels']:
      self.data['mass_per_crystal'] = self.data['mass_per_crystal'] + self.data['scf_settings'][l][1]    
      
    # CELL VOLUMES  
    for i in range(len(self.data['crystals'])):
      self.data['crystals'][i]['cell_volume'] = pwscf_output.cell_volume(self.data['crystals'][i]['alat'], self.data['crystals'][i]['cell_parameters'][:,:])  
      self.data['crystals'][i]['cell_density'] = pwscf_output.cell_density(self.data['crystals'][i]['cell_volume'], self.data['mass_per_crystal'])
      

    

  def next_line(self, n, data):
    if(n < len(data)):
      line = data[n].strip()
      line_uc = line.upper()
      n = n + 1
      return n, line, line_uc
    else:
      n = n + 1
      return n, None, None
    
        
  def store(self, store, line, field, n=0):
    l, f = pwscf_output.read_line(line, field)  
    if(l != False):
      self.data[store] = f[n]

  #  Run as it's own program
  def run(self):
    self.reset()

    option = ""
    file_name = ""

    if(len(sys.argv) > 1 and sys.argv[1] is not None):
      option = sys.argv[1]

    if(len(sys.argv) > 2 and sys.argv[2] is not None):
      file_name = sys.argv[2]

    if(option.lower().strip() == "" or option.lower().strip() == "interactive"):
      self.menu()
      exit()
    elif(option.lower().strip() == "quiet"):
      print("Quiet")
    else:
      return 0






#################################
# READ/LOAD input file
#################################

  def load_from_file(self, file_name):
    # Init variable
    file_data = ""

    # Read it in line by line
    fh = open(file_name, "r")
    for file_row in fh:
      file_data = file_data + file_row.strip() + '\n'

    return file_data




#################################
# Get
#################################

  def get_nat(self): 
    return self.atom_count

  def get_alat(self):
    return self.data['alat']
    
  def get_volume(self):
    return self.data['scf_settings']['volume'] 
    
  def get_volume_per_atom(self):
    return self.data['scf_settings']['volume'] / self.atom_count
  
  def get_total_energy(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy']  
    else:
      return self.data['results'][n]['energy']  
    
  def get_energy_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['energy'] / self.atom_count  
    else:
      return self.data['results'][n]['energy'] / self.atom_count 
    
  def get_total_force(self, n = None):
    if(n == None):
      return self.data['results'][-1]['total_force']  
    else:
      return self.data['results'][n]['total_force']
    
  def get_force_per_atom(self, n=None):
    if(n == None):
      return self.data['results'][-1]['total_force'] / self.atom_count  
    else:
      return self.data['results'][n]['total_force'] / self.atom_count  
  
  def get_density(self):
    return self.data['density']  
    
  def get_cell_parameters(self):
    cp = ['alat', 
          [str(self.data['crystal_calc'][0,0]), str(self.data['crystal_calc'][0,1]), str(self.data['crystal_calc'][0,2])], 
          [str(self.data['crystal_calc'][1,0]), str(self.data['crystal_calc'][1,1]), str(self.data['crystal_calc'][1,2])], 
          [str(self.data['crystal_calc'][2,0]), str(self.data['crystal_calc'][2,1]), str(self.data['crystal_calc'][2,2])]]
    return cp

  # Return relaxed unit vector
  def get_cell_array(self):
    return self.data['crystal_calc']

  # return alat and normalised unit vector
  def get_norm_relaxed(self):
    alat = self.data['crystals'][-2]['alat_adj']
    cp = self.data['crystals'][-2]['cell_parameters_adj']
    return alat, cp

  # Get stress
  def get_stress(self):
    return self.data['stress']
    
  def get_stress_sum(self):
    return self.data['stress_sum']
  

  def get_job_done(self):
    return self.data['job_done']
    
  def get_job_error(self):
    return self.data['error']
    
  def get_job_converged(self):
    return self.data['converged']

  def get_ok(self):
    return self.data['ok']

  # GET - VC-RELAXED

  def get_vc_relax_alat(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['alat_adj']
    else:
      return None

  def get_vc_relax_cp(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_parameters_adj']
    else:
      return None

  def get_vc_relax_volume(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_volume']
    else:
      return None

  def get_vc_relax_density(self):
    if(self.data['type']=="VC-RELAX"):
      return self.data['crystals'][-1]['cell_density']
    else:
      return None
      
  def get_mass_per_crystal(self):
    return self.data['mass_per_crystal']
      
  def get_times(self):
    return self.data['cpu_time'], self.data['wall_time']
      

#################################
# Interactive
#################################

  def menu(self):
    while(True):
      choice = self.print_menu().upper()
      print(choice)
      if(choice == "X"):
        exit()
      elif(choice == "1"):
        self.i_load()
      elif(choice == "2"):
        self.i_display()


  def print_menu(self):
    pwscf_output.header("Menu")
    print("1. Load File")
    print("2. Display File")
    print("X. Exit")
    return input("Choice: ")


  def i_load(self):
    pwscf_output.header("Load Output File")
    file_name = input("Enter file name: ")
    self.load(file_name)
    print("File loaded.")
    input()


  def i_display(self):
    pwscf_output.header("Display File")
    self.output_details()
    input()


  def output_details(self):
    print("Output")
    print("=======================================================================")
    for key in sorted(self.data.keys()):
      value = self.data[key]
      print(key, ":  ", value)
    print("=======================================================================")
    print()
    

  def xyz_evang(self):
    self.xyz_units = 'ev/ang'
    self.stress_units = 'gpa'

  def xyz_stress_gpa(self):
    self.stress_units = 'gpa'

  def make_xyz(self, option=None):
    self.xyz = []
    if(option == None):
      for rn in range(len(self.data['results'])):
        option = rn + 1
    elif(option == -1):
      option = len(self.data['results'])
    else:
      option = (option - 1) % len(self.data['results']) + 1
    self.make_xyz_inner(option)
    return self.xyz
    
    
  def make_xyz_inner(self, option):
    if(len(self.data['results'])==0):
      return False
    if(len(self.data['crystals'])==0):
      return False  
   
    rn = (option - 1) % len(self.data['results'])
    cn = rn
    if(rn == len(self.data['results']) - 1):
      cn = len(self.data['crystals']) - 1
      
    crystal = self.data['crystals'][cn]
    result = self.data['results'][rn]
    settings = self.data['scf_settings']
    species = settings['atomic_species']
    
     
    # Add new list and set counter n
    self.xyz.append([])
    n = len(self.xyz) - 1
    
    # Add data
    self.xyz[n].append("#ALAT " + str(crystal['alat_adj']))
    self.xyz[n].append("#X " + str(crystal['cell_parameters_adj'][0][0]) + " " + str(crystal['cell_parameters_adj'][0][1]) + " " + str(crystal['cell_parameters_adj'][0][2]))
    self.xyz[n].append("#Y " + str(crystal['cell_parameters_adj'][1][0]) + " " + str(crystal['cell_parameters_adj'][1][1]) + " " + str(crystal['cell_parameters_adj'][1][2]))
    self.xyz[n].append("#Z " + str(crystal['cell_parameters_adj'][2][0]) + " " + str(crystal['cell_parameters_adj'][2][1]) + " " + str(crystal['cell_parameters_adj'][2][2]))
    
    # Just use 1 1 1
    self.xyz[n].append("#C 2 2 2")
    self.xyz[n].append("#RCUT 6.5")
    self.xyz[n].append("#L_UNITS bohr")
    self.xyz[n].append("#E_UNITS ry")
    self.xyz[n].append("#F_UNITS ry/bohr")
    self.xyz[n].append("#S_UNITS kbar")
    
    
    self.xyz[n].append("#E " + str(result['energy']))
     
    if(result['s_on']):
      self.xyz[n].append("#SX " + str(result['stress'][0,0]) + " " +  str(result['stress'][0,1]) + " " + str(result['stress'][0,2]))
      self.xyz[n].append("#SY " + str(result['stress'][1,0]) + " " +  str(result['stress'][1,1]) + " " + str(result['stress'][1,2]))
      self.xyz[n].append("#SZ " + str(result['stress'][2,0]) + " " +  str(result['stress'][2,1]) + " " + str(result['stress'][2,2]))
      
    
    for label in species.keys():
      self.xyz[n].append("#M " + label + " " + str(species[label]['mass']))
    
    for i in range(self.atom_count):
      line = crystal['atomic_labels'][i]
      line = line + " " + str(crystal['crystal_positions'][i,0]) + " " + str(crystal['crystal_positions'][i,1]) + " " + str(crystal['crystal_positions'][i,2])      
      if(result['f_on']):   
        line = line + " " + str(result['forces'][i,0]) + " " + str(result['forces'][i,1]) + " " + str(result['forces'][i,2])
    #"s_on": False,
      self.xyz[n].append(line)
      
      
    # Return
    return self.xyz[n]
    
    
    
    
    
    
    
    
   
  def get_data(self, file=None): 
  
    """
    # Important so store in it's own variable
    self.atom_count = 1    

    # Control
    self.data = {
      "ok": False,
      "job_done": False,
      "error": False,
      "type": None,
      "summary": None,
      "mpi_processes": None,
      "threads_per_mpi_process": None,
      
      "scf_settings": None,      
      "crystals": [],
      "results": [],
      
      "initial_positions": None,   
      "total_energy": None,
      "density_full": None,
      "density": None,
      "stress": numpy.zeros((3,3)),
      "stress_sum": None,      
      "cpu_time": None,
      "wall_time": None,   
      "xyz": [],
    }
    """
    
    out = "##############################################################################################################\n"
    out = out + "atom count:                  " + str(self.atom_count) + "\n"
    out = out + "ok:                          " + str(self.data['ok']) + "\n"
    out = out + "job_done:                    " + str(self.data['job_done']) + "\n"
    out = out + "error:                       " + str(self.data['error']) + "\n"
    out = out + "type:                        " + str(self.data['type']) + "\n"
    out = out + "summary:                     " + str(self.data['summary']) + "\n"
    out = out + "mpi_processes:               " + str(self.data['mpi_processes']) + "\n"
    out = out + "threads_per_mpi_process:     " + str(self.data['threads_per_mpi_process']) + "\n"
    out = out + "scf_settings:                " + str(len(self.data['scf_settings'])) + "\n"
    
    for k in self.data['scf_settings'].keys():    
      out = out + "                             " + k + '  '+ str(self.data['scf_settings'][k]) + "\n"  
    
    out = out + "crystals (count):            " + str(len(self.data['crystals'])) + "\n"
    out = out + "results (count):             " + str(len(self.data['results'])) + "\n"
    i = 0
    for i in range(len(self.data['crystals'])-1):
      out = out + "##############################################################################################################\n"
      c = self.data['crystals'][i]
      n = i + 1
      if(n == 1):
        out = out + "crystal " + str(n) + " (input):\n"         
      elif(n == len(self.data['crystals']) - 1):
        out = out + "crystal " + str(n) + " (relaxed):\n" 
      else:
        out = out + "crystal " + str(n) + ":\n"  
      out = out + "##############################################################################################################\n"
      out = out + "alat:                        " + str(c['alat']) + "\n"  
      out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
      out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
      out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
      out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
      out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 
      if(i<len(self.data['crystals'])-2):
        s = self.data['results'][i]
        out = out + "energy:                      " + str(s['energy']) + "\n"  
        out = out + "total_force:                 " + str(s['total_force']) + "\n"  
        out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
        out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
        out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
        for m in range(self.atom_count):
          out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
      
      out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
      out = out +     "density:                     " + str(c['cell_density']) + "\n" 
       
      
    c = self.data['crystals'][-1]
    s = self.data['results'][-1]
    n = len(self.data['crystals'])
    
    out = out + "##############################################################################################################\n"
    out = out + "crystal " + str(n) + " (final): \n"     
    out = out + "##############################################################################################################\n"
    out = out + "alat:                        " + str(c['alat']) + "\n" 
    out = out + "cp:                          " + str(c['cell_parameters'][0,0]) + " " + str(c['cell_parameters'][0,1]) + " " + str(c['cell_parameters'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters'][1,0]) + " " + str(c['cell_parameters'][1,1]) + " " + str(c['cell_parameters'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters'][2,0]) + " " + str(c['cell_parameters'][2,1]) + " " + str(c['cell_parameters'][2,2]) + "\n" 
    out = out + "alat adjusted:               " + str(c['alat_adj']) + "\n"  
    out = out + "cp adjusted:                 " + str(c['cell_parameters_adj'][0,0]) + " " + str(c['cell_parameters_adj'][0,1]) + " " + str(c['cell_parameters_adj'][0,2]) + "\n" 
    out = out + "                             " + str(c['cell_parameters_adj'][1,0]) + " " + str(c['cell_parameters_adj'][1,1]) + " " + str(c['cell_parameters_adj'][1,2]) + "\n"
    out = out + "                             " + str(c['cell_parameters_adj'][2,0]) + " " + str(c['cell_parameters_adj'][2,1]) + " " + str(c['cell_parameters_adj'][2,2]) + "\n" 

    out = out + "energy:                      " + str(s['energy']) + "\n"  
    out = out + "total_force:                 " + str(s['total_force']) + "\n"  
    out = out + "stress:                      " + str(s['stress'][0,0]) + " " + str(s['stress'][0,1]) + " " + str(s['stress'][0,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][1,0]) + " " + str(s['stress'][1,1]) + " " + str(s['stress'][1,2]) + "\n"  
    out = out + "stress:                      " + str(s['stress'][2,0]) + " " + str(s['stress'][2,1]) + " " + str(s['stress'][2,2]) + "\n"   
    out = out + "alat positions,  crystal positions,  forces:   " + "\n" 
    for m in range(self.atom_count):
      out = out + "                             " + str(c['atomic_positions'][m,0]) + " " + str(c['atomic_positions'][m,1]) + " " + str(c['atomic_positions'][m,2]) + "    " + str(c['crystal_positions'][m,0]) + " " + str(c['crystal_positions'][m,1]) + " " + str(c['crystal_positions'][m,2]) + "    " + str(s['forces'][m,0]) + " " + str(s['forces'][m,1]) + " " + str(s['forces'][m,2]) + "\n"
         
    out = out +     "cell volume:                 " + str(c['cell_volume']) + "\n"  
    out = out +     "density:                     " + str(c['cell_density']) + "\n" 
    
    fh = open(file,'w')
    fh.write(out)
    fh.close()
    
    
    return out
    
        
    

#################################
# Static Methods
#################################


  @staticmethod
  def remove_spaces(input_string):
    return input_string.replace(" ", "")
    
    

  @staticmethod
  def extract(input_string, start=None, end=None, type=None, split=None, trim=False):
    if(start == ""):
      start = None
    if(end == ""):
      end = None
    if(trim):
      input_string = input_string.strip()
    
    # Start/End
    start_n = None
    end_n = None
      
    if(start == None and end == None):   
      start_n = 0
      end_n = len(input_string)
    elif(start == None and end != None):  
      end_l = len(end)   
      start_n = 0
      for n in range(len(input_string)):
        if(input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
    elif(start != None and end == None):  
      start_l = len(start)
      end_n = len(input_string)
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
    else:  
      start_l = len(start)
      end_l = len(end)  
    
      for n in range(len(input_string)):
        if(input_string[n:n+start_l] == start[0:start_l]):
          start_n = n + start_l
        if(start_n != None and input_string[n:n+end_l] == end[0:end_l]):
          end_n = n
          break
        
    # Read
    result = input_string[start_n:end_n].strip()       

      
    # Split
    if(split != None):
      if(split == " "):
        #result = re.sub(r'\s\s+', ' ', result)
        result = pwscf_output.single_spaces(result)
      result = result.split(split)
      for i in range(len(result)):
        if(type.lower() == "f"):
          result[i] = float(result[i])
        elif(type.lower() == "i"):
          result[i] = int(result[i])
        
    else:  
      if(type.lower() == "f"):
        result = float(result)
      elif(type.lower() == "i"):
        result = int(result)
        
    # Return
    return result
      
    
  @staticmethod
  def compare(line, field):
    line = line.strip()
    line = line.upper() 
    
    field = field.strip()
    field = field.upper()
    
    f_len = len(field)
    if(len(line) >= f_len and line[0:f_len] == field[0:f_len]):
      return True
    return False
    
    
  @staticmethod
  def read_line(line, field):
    line = line.strip()
    #line = re.sub(r'\s\s+', ' ', line)    
    #line = re.sub(r'\s=\s', '=', line)
    line = pwscf_output.clean(line)
    line_uc = line.upper() 
    
    field = field.strip()
    #field = re.sub(r'\s\s+', ' ', field)
    #field = re.sub(r'\s=\s', '=', field)
    field = pwscf_output.clean(field)
    field = field.upper()
    
    f_len = len(field)
    if(len(line_uc) >= f_len and line_uc[0:f_len] == field[0:f_len]):
      output = line[f_len:].strip()
      fields = output.split(" ")
      return output, fields      
    return False, False
    
 
  @staticmethod
  def fields(input_string):
    input_string = input_string.strip()
    output_string = ""
    last = None
    for character in input_string:
      if(character != " " or (character == " " and last != " ")):
        output_string += character
    return output_string.split(" ")
    
  @staticmethod
  def check_keyword(line, keyword):
    if(line.upper()[0:len(keyword)] == keyword.upper()):
      return True
    return False


  @staticmethod
  def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')


  @staticmethod
  def header(sub_title=""):
    pwscf_output.clear_screen()
    print("==========================================================")
    print("                    PWscf Input Editor                    ")
    print("==========================================================")
    print()
    print(sub_title)
    print()
    print()
    print()

  @staticmethod
  def process_keyword(str_in):
    str_in = str_in.lower().strip()
    str_in = pwscf_output.remove_spaces(str_in)
    id = None
    keyword = ""
    flag = 0
    for character in str_in:
      if(character == "("):
        id = ""
        flag = 1
      elif(character == ")"):
        flag = 2
      elif(flag == 0):
        keyword += character
      elif(flag == 1):
        id = id + character
    if(id != None):
      try:
        id = int(id)
      except:
        id = None
    return keyword, id  

  @staticmethod
  def add_keyword(keywords, keyword, id, value):
    if(id == None):
      added = False
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          added = True
          keywords[i][1] = keyword
      if(added == False):
        keywords.append([keyword, value])
    else:   
      n = None
      for i in range(len(keywords)):
        if(keywords[i][0] == keyword):
          n = i
          break
      if(n == None):    
        keywords.append([keyword,[None]])
        n = len(keywords) - 1
        
      while(len(keywords[n][1]) < id):
        keywords[n][1].append(None)

      keywords[n][1][id-1] = value  



  @staticmethod
  def make_line(key, value):
    output = ""
    if(value != None):
       if(isinstance(value, (list,))):
         for i in range(len(value)):
           if(value[i] != None):
             output += key + "(" + str(i+1) + ") = " + value[i] + ", \n"                
       else:
         output += key + " = " + value + ", \n"   
    return output    

  @staticmethod
  def coord_format(float_in):
    pad = "              "
    value = str(round(float_in, 6)).strip()
    return value
    
  @staticmethod
  def label_format(label):  
    pad = "              "
    label = label.strip()
    return label
    
    

  @staticmethod
  def is_zero(arr):
    for i in range(arr.shape[0]):
      for j in range(arr.shape[1]):
        if(arr[i, j] != 0.0):
          return False
    return True
    
    
  @staticmethod
  def clean(str_in):  
    str_out = ""
    l = len(str_in)
    for i in range(l):
      # Last, Next, This
      if(i == 0):
        last = None
      else:
        last = str_in[i-1]
      if(i < (l-1)):
        next = str_in[i+1]
      else:  
        next = None
      char = str_in[i]
    
      # Check
      ok = True
      if(last == " " and char == " "):
        ok = False
      elif(last == "\n" and char == "\n"):
        ok = False
      elif(last == "\n" and char == " "):
        ok = False
      elif(char == " " and next == "\n"):
        ok = False
      elif(last == "=" and char == " "):
        ok = False
      elif(char == " " and next == "="):
        ok = False
        
      # Add to string  
      if(ok):
        str_out += char
    return str_out    
    
    
    
  @staticmethod
  def electron_string(str_in):
    arr = str_in.split("(up:")
    e = arr[0]
    if(len(arr) == 1):
      return e.strip(), None, None
    if(len(arr)==2):
      arr_b = arr[1].split(", down:")
      eu = arr_b[0]
      arr_c = arr_b[1].split(")")
      ed = arr_c[0]
      return e.strip(), eu.strip(), ed.strip()
  
    print("TEST")
    return "","",""
    
    
    
  @staticmethod
  def single_spaces(str_in):
    str_out = ""
    last = None
    for char in str_in:
      if(char != " " or (char == " " and last != " ")):
        str_out = str_out + char
      last = char
    return str_out
    
    
  @staticmethod
  def cell_volume(alat, cp):
    cp_alat = numpy.zeros((3,3,),)
    cp_alat[:,:] = alat * cp[:,:]
    v = numpy.dot(cp_alat[0,:],numpy.cross(cp_alat[1,:], cp_alat[2,:]))
    return v
    
    
  @staticmethod
  def cell_density(v, mpc):
    v_m3 = v * 1.48036e-31
    m = mpc * 1.66054E-027
    rho = m / v_m3
    return rho
    
    
     
"""
        
  def aaa():
    
    # Load
    ###################################
    n = 0
    counter = 0
    while(n < len(data)):
      n, line, line_uc = self.next_line(n, data)
      if(line != ""):
        counter += 1
        if(counter == 1):
          self.data['summary'] = line
        else:
          if(pwscf_output.compare(line, "Number of MPI processes:")):
            self.data['mpi_processes'] = pwscf_output.extract(line, ":", "", "i") 
        
          if(pwscf_output.compare(line, "bravais-lattice index     =")):
            self.data['bravais_lattice_index'] = pwscf_output.extract(line, "=", "", "i") 
            
          if(pwscf_output.compare(line, "lattice parameter (alat)  =")):
            self.data['alat'] = pwscf_output.extract(line, "=", "a.u.", "f")  
            
          if(pwscf_output.compare(line, "unit-cell volume          =")):
            self.data['volume'] = pwscf_output.extract(line, "=", "(a.u.)^3", "f")     
            
          if(pwscf_output.compare(line, "number of atoms/cell      =")):
            self.data['nat'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of atomic types    =")):
            self.data['types'] = pwscf_output.extract(line, "=", "", "i")  
            
          if(pwscf_output.compare(line, "number of electrons       =")):
            str_e = pwscf_output.extract(line, "=", "", "s")
            e, eu, ed = pwscf_output.electron_string(str_e)
            self.data['electrons'] = e
            self.data['electrons_up'] = eu
            self.data['electrons_down'] = ed
          
            
          if(pwscf_output.compare(line, "number of Kohn-Sham states=")):
            self.data['ks_states'] = pwscf_output.extract(line, "=", "", "i")   
            
          if(pwscf_output.compare(line, "kinetic-energy cutoff     =")):
            self.data['ecutwfc'] = pwscf_output.extract(line, "=", "Ry", "f")  
            
          if(pwscf_output.compare(line, "charge density cutoff     =")):
            self.data['ecutrho'] = pwscf_output.extract(line, "=", "Ry", "f")   
        
          if(pwscf_output.compare(line.strip(), "crystal axes:") and pwscf_output.is_zero(self.data['crystal_in'])):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_in'][j,:] = fields  
              self.data['crystal_calc'][j,:] = fields  
          
          if(pwscf_output.compare(line.strip(), "crystal axes:")):            
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)
              fields = pwscf_output.extract(line, "= (", ")", "s", " ")
              self.data['crystal_calc'][j,:] = fields            
        
          
        
        
          if(pwscf_output.compare(line, "!    total energy")):
            self.data['total_energy'] = pwscf_output.extract(line, "=", "Ry", "f")
            
          if(pwscf_output.compare(line, "Total force =")):
            self.data['total_force'] = pwscf_output.extract(line, "=", "T", "f")
            
            
            
          if(pwscf_output.compare(line, "total   stress  (Ry/bohr**3)")):        
            self.data['stress_sum'] = 0.0
            for j in range(3):              
              n, line, line_uc = self.next_line(n, data)   
              fields = pwscf_output.extract(line, "", "", "f", " ", True)  
              self.data['stress'][j,0] = fields[3] 
              self.data['stress'][j,1] = fields[4] 
              self.data['stress'][j,2] = fields[5]
              self.data['stress_sum'] = self.data['stress_sum'] + abs(fields[0]) + abs(fields[1]) + abs(fields[2])
            
            
#                  "stress": numpy.zeros((3,3)),
#      "stress_sum": None,
            
          if(pwscf_output.compare(line, "density = ")):
            self.data['density_full'] = pwscf_output.extract(line, "=", "", "s")
            self.data['density'] = pwscf_output.extract(line, "=", "g/cm^3", "f")
          
          
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['cpu_time'] = pwscf_output.extract(line, ":", "CPU", "s")
            
          if(pwscf_output.compare(line, "PWSCF        :")):
            self.data['wall_time'] = pwscf_output.extract(line, "CPU", "WALL", "s")
          
          if(pwscf_output.compare(line, "JOB DONE.")):
            self.data['job_done'] = True
            
                    
          if(pwscf_output.compare(line, "Exit code:")):
            self.data['error'] = True  
"""