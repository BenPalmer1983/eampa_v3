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


def main():
  
  # RECORD START TIME
  g.times['start'] = time.time()
  
  # SET UP DIRS
  now = datetime.datetime.now()
  date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
  print(date_time)	
  wd_prefix = now.strftime('wd/%Y%m%d_%H%M%S')
  
  # Set wd
  g.dirs['wd'] = wd_prefix
  
  # Set sub dirs
  if('input' not in g.sub_dirs.keys()):
    g.sub_dirs['input'] = 'input'  
  for sd in g.sub_dirs.keys():
    g.dirs[sd] = g.dirs['wd'] + '/' + g.sub_dirs[sd]
  
  # MAKE DIRS
  for d in g.dirs.keys():
    dir = g.dirs[d]
    std.make_dir(dir)
  
  # OPEN LOG
  g.log_fh = open(g.dirs['log'] + '/main.log', 'w')
  g.log_fh.write('###########################################################################\n')
  g.log_fh.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n')
  g.log_fh.write('###########################################################################\n')
  g.log_fh.write('\n')
  g.log_fh.write('Script: ' + str(sys.argv[0]) + '\n')
  if(len(sys.argv)>1):
    run_program = False
    try:
      g.inp = read_config.read_file(sys.argv[1])
      g.log_fh.write('Loaded: ' + str(sys.argv[1]) + '\n')
      run_program = True
      
      # Copy input
      std.copy(sys.argv[0], g.dirs['input'])
      std.copy(sys.argv[1], g.dirs['input'])
    except:
      g.log_fh.write('Unable to load, exiting\n')
      
    # RUN
    if(run_program):
      eampa.run()
    

  # CLOSE LOG
  g.times['end'] = time.time()
  g.log_fh.write('\n')
  g.log_fh.write('###########################################################################\n')
  g.log_fh.write('Duration: ' + str(g.times['end'] - g.times['start']) + '\n')
  g.log_fh.write('###########################################################################\n')
  g.log_fh.close()



# Run
main()
