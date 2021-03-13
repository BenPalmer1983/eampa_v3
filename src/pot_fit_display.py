

class pf_display:

  display_header_set = False
  last_stage = None

  def clear():    
    os.system('cls' if os.name == 'nt' else 'clear') 
  
  def print_line(w=140):
    for i in range(w):
      print("#", end="")
    print()

  def output():
    try:
      output = int(g.inp['display']['output'])
    except:
      output = 1    
    if(output == 1):
      pf_display_simple.output()
    #if(output == 2):
    #  display.output_2()
    #if(output == 3):
    #  display.output_3()
    #  #print(g.rss)



      
 
 
  @staticmethod
  def pad_r(inp, p=12, r=None):
    if(inp == None):
      return ""  
    if(r is not None):  
      try:
        inp = numpy.round(inp, r)
      except:
        pass
    out = str(inp).strip()  
    while(len(out)<p):
      out = out + " "      
    return out[0:p]

  @staticmethod
  def pad_l(inp, p=12, r=None):
    if(inp == None):
      return ""  
    if(r is not None):  
      try:
        inp = numpy.round(inp, r)
      except:
        pass  
    out = str(inp).strip()  
    while(len(out)<p):
      out = " " + out     
    return out[0:p]
    
  @staticmethod
  def bar(p, w=25):
    out = ''
    p_num = p
    p = p * (25 / 100)
    for i in range(w):
      if(i <= p):
        out = out + '#'
      else:
        out = out + '_'
    out = out + '   ' + str(p_num) + '%'
    return out
    