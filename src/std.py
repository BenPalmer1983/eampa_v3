
import re
import os
import numpy
import shutil

class std:

  @staticmethod
  def file_to_list(file_name, clean=False):
    # Init variable
    file_data = []
    # Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(clean):
        line = line.strip()
        if(line != ""):
          file_data.append(line)          
      else:
        file_data.append(line[0:-1])
    # Return
    return file_data
    
  @staticmethod
  def split_fields(line, sep=" "):
    out = line.split(sep)
    key = out[0]
    value = out[1]
    value_out = ''    
    indata = False
    for char in value:
      if(indata and char != '"'):
        value_out = value_out + char
      elif(indata and char == '"'):
        indata = False
      elif(not indata and char == '"'):
        indata = True
    return key, value_out
    
  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out   
    
  @staticmethod
  def to_fields(line, sep=" "):
    out = []
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == sep and last_char == sep)):
        if(char == sep):
          temp = temp.strip()
          if(temp != ""):
            out.append(temp)
            temp = ''
        else:
          temp = temp + char
    
    temp = temp.strip()
    if(temp != ""):
      out.append(temp)      
    return out    
    
    
  @staticmethod
  def make_dir(dir):
    dirs = dir.split("/")
    try:
      dir = ''
      for i in range(len(dirs)):
        dir = dir + dirs[i]
        if(not os.path.exists(dir) and dir.strip() != ''):
          os.mkdir(dir) 
        dir = dir + '/'
      return True
    except:
      return False
      
      
      
  @staticmethod
  def copy(src, dest):  
    try:
      std.make_dir(dest)
      src_f = src.split("/")
      if(os.path.isdir(src)):
        shutil.copytree(src, dest + '/' + src_f[-1])      
      else:
        shutil.copyfile(src, dest + '/' + src_f[-1])
      return True
    except:
      return False
      
  @staticmethod
  def path(dir, file):  
    dir = dir.strip()
    file = file.strip()
    if(dir == ''):
      return file
    else:
      if(dir[-1] == '/'):
        return dir + file
      else:      
        return dir + '/' + file
    
  @staticmethod
  def remove_comments(content):
    data = ''
    i = 0
    for line in content:
      if(i > 0):
        data += '\n'
      data += line
      i = i + 1
    out = ''
    indata = 0
    incomment = 0
    for i in range(len(data)):
      # Get char and next char
      char = data[i]
      next = None
      prev = None
      if(i < len(data)-1):
        next = data[i + 1]
      if(i > 0):
        prev = data[i - 1]
      # If in '  '  
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      # If in "  "  
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0):
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 1 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "!"):
          incomment = 2
        elif(incomment == 2 and char == "\n"):
          incomment = 0
        if(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
    return out.split("\n")    
    
    
  # Remove comments from a block of data/text  
  @staticmethod
  def remove_comments_data(data):
    out = ""
    n = 0
    inquotes = 0
    incomment = 0
    while n < len(data):
      # Get char and next char
      char = data[n]
      next = None
      prev = None
      if(n < len(data)-1):
        next = data[n + 1]
      if(n > 0):
        prev = data[n - 1]
        
      # If in '  '  
      if(inquotes == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(inquotes == 1 and char == "'" and last_char != "\\"):
        out = out + char
        inquotes = 0
      # If in "  "  
      elif(inquotes == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(inquotes == 2 and char == '"' and last_char != "\\"):
        out = out + char
        inquotes = 0
      # If not inside quotes  
      elif(inquotes == 0):
        # Comment on a line
        if(incomment == 0 and char == "/" and next == "/"):
          incomment = 1
        elif(incomment == 0 and char == "!"):
          incomment = 1
        elif(incomment == 0 and char == "#"):
          incomment = 1    
        # Comment on line close
        elif(incomment == 1 and char == "\n"):
          incomment = 0
        # Comment block
        elif(incomment == 0 and char == "/" and next == "*"):
          incomment = 3
        elif(incomment == 3 and prev == "*" and char == "/"):
          incomment = 0
        elif(incomment == 0):
          out = out + char  
      # Increment counter
      n = n + 1
    return out        

      
  # Single spaces, tabs to spaces        
  @staticmethod
  def prep_data(content):
    out = []
    for line in content:
      line_new = std.prep_data_line(line)
      if(line_new != ''):
        out.append(line_new)
    return out  
      
  @staticmethod
  def prep_data_line(line): 
    temp = ''
    indata = 0
    last_char = None
    for char in line:
      if(char == '\t'):
        char = ' '
      if(indata == 1 and char != "'" and last_char != "\\"):
        temp = temp + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        temp = temp + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        temp = temp + char
        indata = 0
      elif(indata == 0 and not (char == ' ' and last_char == ' ')):
        temp = temp + char       
      last_char = char  
    return temp.strip()    
    
    
  @staticmethod
  def remove_quotes(inp): 
    if(isinstance(inp, list)):    
      for i in range(len(inp)):
        inp[i] = std.remove_quotes(inp[i])        
      return inp
    else:
      inp = inp.strip()
      if(inp[0] == '"' and inp[-1] == '"'):
        return inp[1:-1]
      if(inp[0] == "'" and inp[-1] == "'"):
        return inp[1:-1]
      return inp
      
    
    
    
  @staticmethod
  def config_file_to_list(file_name):
    # Init variable
    file_data = []
    # Read it in line by line
    fh = open(file_name, "r")
    for line in fh:
      if(line.strip() != ""):
        line = line.strip()
        line = std.remove_comments(line)
        line = std.prep_data_line(line)
        fields = std.to_fields(line)
        file_data.append(fields)         
    # Return
    file_data = std.remove_quotes(file_data)
    return file_data
    
    
  @staticmethod
  def get_dir(file_path):
    directory = ''
    read = False
    for i in range(len(file_path)):
      if(read):
        directory = file_path[-1-i] + directory
      if(file_path[-1-i] == "/"):
        read = True
    return directory
  
    
  @staticmethod
  def read_csv(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
      # Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line
      fh.close()
      # Remove comments
      file_data = std.remove_comments_data(file_data)
      # Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
    return data
     
    
  @staticmethod
  def read_csv_array(filename, sep=","):
    data = []
    if(os.path.isfile(filename)):
      # Read from file into memory
      fh = open(filename, 'r')
      file_data = ""
      for line in fh:
        file_data = file_data + line.strip() + '\n'
      fh.close()
      # Remove double spaces
      file_data = std.one_space(file_data)
      # Remove comments
      file_data = std.remove_comments_data(file_data)
      # Read Data
      lines = file_data.split("\n")
      for line in lines:
        line = line.strip()
        if(line != ""):
          data.append(line.split(sep))  
          
      lst = []
      for i in range(len(data)):
        row = []
        for j in range(len(data[0])):
          try:
            row.append(float(data[i][j]))
          except:
            pass
        lst.append(row)
        
          
      arr = numpy.zeros((len(lst),len(lst[0]),),)
      for i in range(len(lst)):
        for j in range(len(lst[0])):
          arr[i,j] = float(lst[i][j])
      return arr
    return None
    
  @staticmethod
  def write_csv(filename, arr):  
    fh = open(filename, 'w')
    for i in range(len(arr)):
      for j in range(len(arr[i])):
        fh.write(std.float_padded(arr[i,j],8) + " ")
      fh.write("\n")
    fh.close()
  
    
    
  @staticmethod
  def option(input):
    input = input.strip().upper()
    if(input[0:1] == "Y"):
      return True
    elif(input[0:2] == "ON"):
      return True
    elif(input[0:1] == "T"):
      return True
    else:
      return False
    
    
  @staticmethod
  def float_padded(inp, pad=7):
    out = float(inp)
    out = round(out, pad-3)
    out = str(out)  
    while(len(out)<pad):
      out = out + " "      
    return out[0:pad]
    
  @staticmethod
  def write_file_line(fh, title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    fh.write(line)  
  
  
  @staticmethod
  def print_file_line(title, title_pad, fields, field_pad):
    if(type(fields) == numpy.ndarray):
      t = fields
      fields = []
      for ti in t:
        fields.append(ti)    
    elif(type(fields) != list):
      fields = [fields]
    
    line = str(title)
    while(len(line)<title_pad):
      line = line + ' '
    for f in fields:
      f_str = str(f)
      while(len(f_str)<field_pad):
        f_str = f_str + ' '
      line = line + f_str + ' '
    line = line + '\n'
    print(line,end='')  
    
    
  def pad(inp, width):
    out = str(inp)
    while(len(out)<width):
      out = out + " "
    return out
    
  def mem_value(strin):
    num = '0123456789.'
    strin = strin.upper().strip()
    value = ''
    unit = ''
    for c in strin:
      if(c not in num):
        break
      else:
        value = value + c
    for c in strin:
      if(c not in num):
        unit = unit + c
     
    value = float(value)
    if(unit == 'B'):
      return value
    if(unit == 'KB'):
      return 1000 * value
    if(unit == 'MB'):
      return 1000000 * value
    if(unit == 'GB'):
      return 1000000000 * value
    
    
    
    
  @staticmethod
  def dict_to_str(d, pre=''):
    out = ''
    for k in sorted(d):
      if(type(k) == dict or type(k) == list):
        out = out + '\n' + std.dict_to_str(k, pre + '  ')      
      else:
        out = out + pre + str(k) + ': ' + str(d[k]) + '\n'
    return out
    
    
    
    