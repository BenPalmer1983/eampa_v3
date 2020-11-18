################################################################
#    Main Program
#
#
#
#
################################################################

import os
import time
import datetime
import re
import sys
import shutil
from globals import g
from std import std
from read_config import read_config
from eampa import eampa

class main():

  log_info_count = 0

  def start():
  
    # RECORD START TIME
    g.times['start'] = time.time()
  
    now = datetime.datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print(date_time)	
  
    # OPEN LOG
    g.log_fh = open('job.log', 'w')
    main.log_hr()
    main.log(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
    main.log_hr()
    main.log_br()
    main.log_info()
    main.log('Script: ' + str(sys.argv[0]))
    if(len(sys.argv)>1):
      run_program = False
      try:
        g.inp = read_config.read_file(sys.argv[1])
        main.log('Loaded: ' + str(sys.argv[1]))
        run_program = True
      
        # Copy input
        std.copy(sys.argv[0], g.dirs['input'])
        std.copy(sys.argv[1], g.dirs['input'])
      except:
        main.log('Unable to load, exiting\n')
      
      # RUN
      if(run_program):
        eampa.run()
        
      # End and close log
      main.end()
      
      
  def log(line='', end='\n'): 
    main.log_info_count = main.log_info_count + 1
    if(main.log_info_count == 100):
      main.log_info_count = 0
      main.log_info()
  
    lines=line.split("\n")
    for line in lines:
      g.log_fh.write(str("{:.3E}".format(time.time() - g.times['start'])) + ' ###   ' + line + end) 
      
  def log_hr(): 
    g.log_fh.write('#############################################################################\n') 
    
  def log_br(): 
    g.log_fh.write('\n') 
        
  def log_info(): 
    g.log_fh.write('####TIME#############LOG#####################################################\n') 
    
  def log_title(title=''): 
    g.log_fh.write('#############################################\n') 
    titles = title.split("\n")
    for title in titles:
      g.log_fh.write('  ' + title + '\n') 
    g.log_fh.write('  Time: ' + str("{:.3E}".format(time.time() - g.times['start'])) + '\n') 
    g.log_fh.write('#############################################\n') 
  
   
  def end(): 

    # CLOSE LOG
    g.times['end'] = time.time()
    main.log_br()
    main.log_hr()
    main.log('Duration: ' + str(g.times['end'] - g.times['start']))
    main.log_hr()
    g.log_fh.close()
    exit()



# Run
main.start()