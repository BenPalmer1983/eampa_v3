######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################

import numpy



class relax_calc:

  def run():  
  
    print("Relax") 
    
    relax.init()
    potential.relax_add_potentials()
    
    c = g.configs['configs'][0]  
    #relax_calc.randomise_coords(c)
    relax.add_config( 
                     7.0,
                     c['alat'], 
                     c['uv'], 
                     c['coords_label_id'], 
                     c['coords']
                    )  
               
         
    relax.run(25, 0.07, -0.9)  
               
    relax_calc.save_xyz()
    
    
      
  def md():  
    relax.run_md(1000, 0.01)                
    relax_calc.save_xyz_md()
    
  def randomise_coords(c):     
    for ci in range(len(c['coords'])):
      r = numpy.random.rand(3)
      c['coords'][ci, :] = c['coords'][ci, :]  + (0.5- r) * 0.04
    
      
  def save_xyz():
    dir = g.dirs['wd'] + "/relax"
    std.make_dir(dir)
    fh = open(dir + '/relaxed.xyz', 'w')
    fh.write(str(relax.atom_count) + "\n")
    fh.write("\n")
    for i in range(relax.atom_count):
      fh.write(str(labels.get(relax.labels[i])) + " ")
      fh.write(str(relax.coords[i, 3]) + " ")
      fh.write(str(relax.coords[i, 4]) + " ")
      fh.write(str(relax.coords[i, 5]) + " ")
      fh.write("\n")      
    fh.close() 
    fh = open(dir + '/relaxed_crystal.xyz', 'w')
    fh.write(str(relax.atom_count) + "\n")
    fh.write("\n")
    for i in range(relax.atom_count):
      fh.write(str(labels.get(relax.labels[i])) + " ")
      fh.write(str(relax.coords[i, 0]) + " ")
      fh.write(str(relax.coords[i, 1]) + " ")
      fh.write(str(relax.coords[i, 2]) + " ")
      fh.write("\n")      
    fh.close() 
  
  def save_xyz_md():
    dir = g.dirs['wd'] + "/md"
    std.make_dir(dir)
    
    fh = open(dir + '/md.xyz', 'w')
    for n in range(relax.md_steps + 1):
      fh.write(str(relax.atom_count) + "\n")
      fh.write(str(n) + "\n")
      for i in range(relax.atom_count):
        fh.write(str(labels.get(relax.labels[i])) + " ")
        fh.write(str(relax.md_xyz[n, i, 0]) + " ")
        fh.write(str(relax.md_xyz[n, i, 1]) + " ")
        fh.write(str(relax.md_xyz[n, i, 2]) + " ")
        fh.write("\n")
      
    fh.close()           
               
    #relax.make_nl()                
    #relax.force_calc()        
    #print(relax.nl_count)
    #print(relax.config_forces[0:relax.atom_count, :])
    
    # Setup ES
    #es.init()
    #potential.es_add_potentials()
    #es_calc.es_add()