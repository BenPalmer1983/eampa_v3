import sys
import numpy


class config:

  @staticmethod
  def run():
    print("Expand Config")
    
    if(len(sys.argv) != 3):
      print("python3 config.py input.conf output.conf")
      exit()

    infile = str(sys.argv[1])
    outfile = str(sys.argv[2])

    fh = open(infile, 'r')
    c = []
    infile_lines = []
    for line in fh:
      line = line.strip()
      if(line != ""):
        infile_lines.append(line)
      if(line[0:2].upper() == '#C'):
        f = config.one_space(line).split(" ")
        cx = int(f[1])
        cy = int(f[2])
        cz = int(f[3])
      elif(line[0:2].upper() == '#X'):
        x = config.one_space(line).split(" ")
      elif(line[0:2].upper() == '#Y'):
        y = config.one_space(line).split(" ")
      elif(line[0:2].upper() == '#Z'):
        z = config.one_space(line).split(" ")
      elif(line[0] != "#"):
        f = config.one_space(line).split(" ")
        if(len(f) == 4 or len(f) == 7):
          c.append(f)
    fh.close()

    uv = numpy.asarray([[float(x[1]),float(x[2]),float(x[3])],[float(y[1]),float(y[2]),float(y[3])],[float(z[1]),float(z[2]),float(z[3])]])
    cm = numpy.asarray([[float(cx),0.0,0.0],[0.0,float(cy),0.0],[0.0,0.0,float(cz)]])
    uv_new = numpy.matmul(cm, uv)
    
    print(infile_lines)
    fh = open(outfile, 'w')
    for line in infile_lines:
      if(line[0:2].upper() == '#C'):
        fh.write("#C 1 1 1\n")
      elif(line[0:2].upper() == '#X'):
        fh.write("#X " + str(uv_new[0,0]) + " " + str(uv_new[0,1]) + " " + str(uv_new[0,2]) + "\n")
      elif(line[0:2].upper() == '#Y'):
        fh.write("#Y " + str(uv_new[1,0]) + " " + str(uv_new[1,1]) + " " + str(uv_new[1,2]) + "\n")
      elif(line[0:2].upper() == '#Z'):
        fh.write("#Z " + str(uv_new[2,0]) + " " + str(uv_new[2,1]) + " " + str(uv_new[2,2]) + "\n")
      elif(line[0:1].upper() == '#'):
        fh.write(line + "\n")

    for i in range(cx):
      for j in range(cy):
        for k in range(cz):
          for cn in c: 
            c_x = (i + float(cn[1])) / cx
            c_y = (j + float(cn[2])) / cy
            c_z = (k + float(cn[3])) / cz
            fh.write(cn[0] + " " + str('{:16.5e}'.format(c_x)) + " " + str('{:16.5e}'.format(c_y)) +  " " + str('{:16.5e}'.format(c_z)))
            if(len(cn) == 7):
              fh.write(" " + str('{:16.5e}'.format(float(cn[4]))) + " " + str('{:16.5e}'.format(float(cn[5]))) +  " " + str('{:16.5e}'.format(float(cn[6]))))
            fh.write('\n')

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


config.run()

