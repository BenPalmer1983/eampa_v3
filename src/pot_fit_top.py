class pf_top:

  def save(dir, filename):

    print("Saving Top Parameters")
    print(dir + "/" + filename)
    std.make_dir(dir)
    fh = open(dir + "/" + filename, 'w')
    for i in range(len(g.top_parameters)):
      i_str = str(i+1)
      while(len(i_str)<8):
        i_str = "0" + i_str
      fh.write("---START---" + "\n")
      fh.write(i_str + "\n")
      fh.write("PARAMETERS \n")
      fh.write(str(len(g.top_parameters[i][1])) + "\n")
      for j in range(len(g.top_parameters[i][1])):
        fh.write(str(g.top_parameters[i][1][j]) + "\n")
      fh.write("RSS " + str(g.top_parameters[i][0]) + '\n')
      fh.write("BP MaxDensity " + str(g.top_parameters[i][4]) + '\n')
      fh.write("EFS MaxDensity " + str(g.top_parameters[i][5]) + '\n')
      try:
        bp_calculations = g.top_parameters[i][3]['bp_calculations']
        for bi in range(len(bp_calculations)):
          fh.write("a0 " + str(bp_calculations[bi]['a0']) + "\n")
          fh.write("e0 " + str(bp_calculations[bi]['e0']) + "\n")
          fh.write("b0 " + str(bp_calculations[bi]['b0']) + "\n")
      except:
        pass
      fh.write("---END---" + "\n")

    fh.close()


  def load(dir, filename):
    if(int(g.fitting['load_top_parameters']) == 0):
      return []
    print(g.fitting['load_top_parameters'])
    top_p = []
    fh = open(dir + "/" + filename, 'r')
    
    f = []
    for line in fh:
      f.append(line.strip())
    fh.close()

    i = 0
    while(i<len(f)):
      if(f[i] == "---START---"):
        i = i + 3
        pn = int(f[i])
        p = numpy.zeros((pn,),)
        for n in range(pn):
          i = i + 1
          p[n] = float(f[i])
      elif(f[i] == "---END---"):
        top_p.append(p)
      i = i + 1
      if(len(top_p) == g.fitting['load_top_parameters']):
        break

    return top_p 
