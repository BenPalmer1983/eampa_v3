#!/bin/python3
################################################################################

#  Combine into one file
#  Changes directory to ../class_files
#  Package main.py
#


import os
import sys
import datetime




class pack:

  def __init__(self):

    # File name
    self.output_file = "main.py"
    if(len(sys.argv) > 1):
      self.output_file = sys.argv[1]
      if(self.output_file[-3:] != ".py"):
        self.output_file = self.output_file + ".py"

    # Record directories
    self.top = os.getcwd()
    os.chdir(self.top + "/src")  
    self.class_files = os.getcwd()

    # Create lists
    self.classes = []
    self.imports = []
    self.file_lines = []
    self.file_lines_out = [] 

    # Run
    self.read_files("main.py")
    self.make()
    self.read_main("main.py")
    self.clean()
    self.output()
    return None


  def read_files(self, file_name):

    # Open to read
    fh = open(file_name, "r")
    for line in fh:

      # Store the IMPORTS
      if(line.strip()[0:6] == "import"):
        line = strf.strip_double_spaces(line.strip())
        imlist = line.split(' ')
        if(len(imlist) >= 2):
          if(line not in self.imports):
            self.imports.append(line)

      # Store the class files
      if(line.strip()[0:4] == "from"):
        line = strf.strip_double_spaces(line.strip())
        imlist = line.split(' ')
        next_file = imlist[1] + ".py"

        if(os.path.isfile(next_file)):
          if(len(imlist) == 4):
            key = imlist[1] + "||" + imlist[3]
            if(key not in self.classes):
              self.classes.append(key)
              self.read_files(imlist[1] + ".py")
        else:
          if(line not in self.imports):
            self.imports.append(line)

      # Store the class files
      if(line.strip()[0:5] == "class"):
        l1 = line.split(":")
        l2 = l1[0].split(" ")       
        key = file_name[0:-3] + "||" + l2[-1]
        if(key not in self.classes):
          self.classes.append(key)
        

    fh.close()



  def make(self):
    horiz = "########################################################################"
    self.file_lines.append("#!/bin/python3")
    self.file_lines.append(horiz)
    for im in self.imports:
      self.file_lines.append(im)
    for key in self.classes:
      fields = key.split("||")
      if(len(fields) == 2):
        self.extract_class(fields[0] + ".py", fields[1].strip())


  def extract_class(self, class_file, class_name):
    class_file = self.class_files + "/" + class_file
    if(os.path.isfile(class_file)):
      fh = open(class_file, "r")
      read = False
      for line in fh:
        if(read == False and line[0:5] == "class"):
          cl = line.strip().split(" ")
          if(len(cl) == 2 and cl[1] == class_name + ":"):
            read = True
            self.file_lines.append("")
            self.file_lines.append("###########################################")
            self.file_lines.append("#  CLASS " + class_name[0:-1])
            self.file_lines.append("###########################################")
        elif(read == True and line[0:5] == "class"):
          break
        if(read):
          self.file_lines.append(line.replace("\n",""))
      fh.close()


  # Finally read main file
  def read_main(self, file_name): 
    before = []
    after = []
    after.append("")
    after.append("###########################################")
    after.append("###########################################")
    after.append("#  MAIN ")
    after.append("###########################################")
    after.append("###########################################")
    after.append("")

    # Open to read
    fh = open(file_name, "r")
    flag = "before"
    for line in fh: 

      if(line.strip()[0:6] == "import" or line.strip()[0:4] == "from"):
        if(flag == "before"):
          flag = "after"
      else:
        if(flag == "before"):
          before.append(line)
        elif(flag == "after"):
          after.append(line)

    #    print(line)  
    #    self.file_lines
    fh.close()

    # Join lists
    self.file_lines = before + self.file_lines + after


  def clean(self):
    for i in range(len(self.file_lines)):
      if(self.file_lines[i].strip()[0:1] == "#"):
        self.file_lines_out.append(self.file_lines[i].strip())
      elif(i == 0 and self.file_lines[i].strip() == ""):
        pass
      elif(i > 0 and self.file_lines[i-1].strip() == "" and self.file_lines[i].strip() == ""):
        pass
      else:
        self.file_lines_out.append(self.file_lines[i])

    
          

  def output(self):
    os.chdir(self.top)
    fh = open(self.output_file, 'w')
    for line in self.file_lines_out:
      fh.write(line + "\n")
    fh.close()


##############################################



class strf:

  @staticmethod
  def pad(data, width=64, align="L"):
    data = str(data)
    in_len = len(data)

    if(in_len<width):
      if(align.upper()[0:1] == "L"):
        for i in range(0, width-in_len):
          data = data + " "

    return data


  @staticmethod
  def strip_double_spaces(string, tabs = True):
    string_out = ""
    last_letter = None
    for letter in string:
      if(tabs and letter == "\t"):
        letter = " "
      if(last_letter != " " or (last_letter == " " and letter != " ")):
        string_out += letter
      last_letter = letter
    return string_out

  @staticmethod
  def parse_line(line):    
    atom_fields = line.split(" ")      
    coords = []
    if(len(atom_fields) == 4):
      # set coords
      for i in range(1,4):
        coords.append(float(atom_fields[i]) )
      return True, atom_fields[0], coords
    return Fale, "", [0,0,0]


  @staticmethod
  def str_compare(str_A, str_B, limit_len=None):  
    # Length of left (str_A)
    if(limit_len is not None and limit_len.upper()[0] == "L"):
      if(len(str_B)<len(str_A)):
        return False
      if(str_A == str_B[0:len(str_A)]): 
        return True
      return False

    # Length of right (str_B)
    if(limit_len is not None and limit_len.upper()[0] == "R"):
      if(len(str_A)<len(str_B)):
        return False
      if(str_A[0:len(str_B)] == str_B): 
        return True
      return False
  
    # No length specified
    if(str_A == str_B): 
      return True
    return False

  @staticmethod
  def str_icompare(str_A, str_B, limit_len=None):  
    str_A = str_A.strip().lower()
    str_B = str_B.strip().lower()
    return strf.str_compare(str_A, str_B, limit_len)






##############################################



# RUN
p = pack()






##############################################


