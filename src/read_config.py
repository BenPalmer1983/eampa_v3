######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

class read_config:
  
  @staticmethod
  def read_file(file_path):
  
    # Input dictionary
    input = {}
  
    # READ DATA
    d = []
    fh = open(file_path, 'r')
    for line in fh:
      line = line.strip()
      if(len(line) > 0 and line[0] != "#"):
        d.append(line)      
    fh.close()
    
    # Count commands
    commands = {}
    for line in d:
      fields = read_config.split_by(line, ' ')
      c = fields[0].lower()
      if(c in commands.keys()):
        commands[c] = commands[c] + 1
      else:
        commands[c] = 1
        
    # Prepare input dictionary    
    for k in commands.keys():
      if(commands[k] == 1):
        input[k] = None
      else:
        input[k] = []
    
    # Read Data into input
    for line in d:
      fields = read_config.split_by(line, ' ')
      fkey = fields[0].lower()
      
      fd_size = {}
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")
        fc = fs[0].lower()
        if(fc in fd_size.keys()):
          fd_size[fc] = fd_size[fc] + 1
        else:
          fd_size[fc] = 1
          
      # Prepare dictionary   
      fd = {} 
      for k in fd_size.keys():
        if(fd_size[k] == 1):
          fd[k] = None
        else:
          fd[k] = []        
        
      for i in range(1, len(fields)):
        f = fields[i]
        fs = f.split("=")     
        fc = fs[0].lower()        
        fs = read_config.split_by(fs[1], ',')         
        fs = read_config.store(fs)
        
        if(fd_size[fc] == 1):
          if(len(fs) == 1):
            fd[fc] = read_config.store(fs[0])
          else:
            fd[fc] = read_config.store(fs)
        else:
          if(len(fs) == 1):
            fd[fc].append(read_config.store(fs[0]))
          else:
            fd[fc].append(read_config.store(fs))
            
      if(commands[fkey] == 1):
        input[fkey] = fd
      else:
        input[fkey].append(fd)  

    return input
        
        
    
  @staticmethod  
  def split_by(line, sep=' ', ignore_double_sep=True):
    last_char = None
    in_quotes = 0
    fields = []
    temp_line = ""
    
    for char in line:
      if(char == "'" and in_quotes == 0 and last_char != "\\"):
        in_quotes = 1
      elif(char == "'" and in_quotes == 1 and last_char != "\\"):
        in_quotes = 0
      elif(char == '"' and in_quotes == 0 and last_char != "\\"):
        in_quotes = 2
      elif(char == '"' and in_quotes == 2 and last_char != "\\"):
        in_quotes = 0
      elif(in_quotes > 0):
        temp_line = temp_line + char
      elif(in_quotes == 0 and char != sep):
        temp_line = temp_line + char
      elif(char == sep and last_char == sep and ignore_double_sep):
        pass
      elif(char == sep):
        fields.append(temp_line)
        temp_line = "" 
    if(temp_line != ""):
      fields.append(temp_line)
    
    return fields
    
    
  @staticmethod
  def store(inp):  
    if(isinstance(inp, list)):
      for i in range(len(inp)):
        try:
          if('.' in inp[i]  or 'e' in inp[i]):
            inp[i] = float(inp[i])
          else:
            inp[i] = int(inp[i])
        except:
          pass
    else:
      try:
        if('.' in inp or 'e' in inp):
          inp = float(inp)
        else:
          inp = int(inp)
      except:
        pass
    return inp
      