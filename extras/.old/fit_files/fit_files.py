import os


class fit_files:


  @staticmethod
  def run():
    file_path = 'files'

    with os.scandir(file_path) as it:
      for entry in it:
        if entry.name.endswith(".pot") and entry.is_file():
          fit_files.make_fit(file_path, str(entry.name))


  @staticmethod
  def make_fit(file_path, file_name):
    
    pot_file_name = file_name.split(".")
    fit_file_name = pot_file_name[0] + ".fit" 

    p = []
    fh = open(file_path + '/' + file_name, 'r')
    for line in fh:
      line = line.strip()
      if(line[0:3] == "#P "):
        line = fit_files.one_space(line)
        r = line.split(' ')
        for i in range(1, len(r), 1):
          p.append(float(r[i]))
    fh.close()


    u = []
    l = []
    for i in range(len(p)):
      u.append(p[i] + 0.1 * p[i])
      l.append(p[i] - 0.1 * p[i])


    fh = open(file_path + '/' + fit_file_name, 'w')
    fh.write("#FIT A \n")
    fh.write("#PS ")
    for pi in p:
      fh.write(str(pi) + " ")
    fh.write("\n")      
    fh.write("#PL ")
    for li in l:
      fh.write(str(li) + " ")
    fh.write("\n")    
    fh.write("#PU ")
    for ui in u:
      fh.write(str(ui) + " ")
    fh.write("\n")   
    fh.close() 



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

fit_files.run()
