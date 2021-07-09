
from globals import g





class setup_dirs:
 
 
  def run():
    main.log_title("Setup Directories")
 
    if(g.wd_type == 1):
      
      # SET UP DIRS
      now = datetime.datetime.now()
      date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
      wd_prefix = now.strftime('wd/%Y%m%d_%H%M%S')
  
      # Set wd
      g.dirs['wd'] = wd_prefix
    
    else:
  
      # Set wd
      g.dirs['wd'] = 'wd/wd'
      
      
    main.log("WD: " + g.dirs['wd'])
    
    
    # Set sub dirs
    if('input' not in g.sub_dirs.keys()):
      g.sub_dirs['input'] = 'input'  
    for sd in g.sub_dirs.keys():
      g.dirs[sd] = g.dirs['wd'] + '/' + g.sub_dirs[sd]
  
    
    
    # For copies of all data to be stored in (potential etc)
    std.make_dir(g.dirs['wd'] + '/data')
    
    
    