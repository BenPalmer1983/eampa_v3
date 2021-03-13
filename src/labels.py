class labels:

  @staticmethod
  def add(label):
    label = label.upper()
    
    # Check for mask
    if(label in g.mask.keys()):
      label = g.mask[label].upper()
    
    if(label in g.labels.keys()):
      return label, g.labels[label]
    else: 
      g.labels[label] = len(g.labels) + 1
      return label, g.labels[label]
     
  @staticmethod
  def output(): 
    if(g.outputs):     
      fh = open(g.dirs['output'] + '/' + 'labels.dat', 'w')   
      for k in g.labels.keys():
        fh.write(str(k) + '  ' + str(g.labels[k]) + '\n')
      fh.close()  
        
  @staticmethod
  def get(l_id):     
    for k in g.labels.keys():
      if(l_id == g.labels[k]):
        return k
    return None
        
       
  @staticmethod
  def add_group(label, group, fn): 
    label = str(label.upper()) 
    group = str(group.upper()) 
    key = label + group
    if(key not in g.grouplabels.keys()):
      g.grouplabels[key] = len(g.grouplabels) + 1
    if(label not in g.groupelement.keys()):
      g.groupelement[label] = [fn]
    else:
      g.groupelement[label].append(fn)
    return key, g.grouplabels[key]
    
  @staticmethod
  def get_group(l_id):     
    for k in g.grouplabels.keys():
      if(l_id == g.grouplabels[k]):
        return k
    return None
        