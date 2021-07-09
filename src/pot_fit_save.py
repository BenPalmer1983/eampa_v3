






class pf_save:
 
  count = 0
  count_name = ''

  def top(name):
    pf_save.count = pf_save.count + 1

    # Get Best
    rss = g.pfdata['rss']['best']
    p = numpy.copy(g.pfdata['p']['best'])
    bp = g.pfdata['bp']['best'].copy()
    bp_known = g.pfdata['bp']['known'].copy()
    efs = g.pfdata['efs']['best'].copy()
    efs_known = g.pfdata['efs']['known'].copy()
 
    # Make Dirs
    c = str(pf_save.count)
    while(len(c)<4):
      c = '0' + c
    name = c + '_' + name
    pf_save.count_name = name

    dir_out = g.dirs['wd'] + '/fitting/saved_potentials/' + name + '/'
    plot_dir = dir_out + 'plots/'
    pot_dir = dir_out + 'pots/'
    std.make_dir(dir_out)
    std.make_dir(plot_dir)
    std.make_dir(pot_dir)


    # Save Top Ten
    fh = open(dir_out + "/top_parameters.txt", 'w')
    for i in range(len(g.top_parameters)):
      i_str = str(i+1)
      while(len(i_str)<8):
        i_str = "0" + i_str
      fh.write('########################' + '\n')
      fh.write('# ' + i_str + '\n')
      fh.write('########################' + '\n')
      fh.write("RSS " + str(g.top_parameters[i][0]) + '\n')
      try:
        bp_calculations = g.top_parameters[i][3]['bp_calculations']
        for bi in range(len(bp_calculations)):
          fh.write("BP " + str(bi) + ":\n")
          fh.write("a0: " + '{:12.4e}'.format(bp_calculations[bi]['a0']) + "   ")
          fh.write("e0: " + '{:12.4e}'.format(bp_calculations[bi]['e0']) + "   ")
          fh.write("b0: " + '{:12.4e}'.format(bp_calculations[bi]['b0_gpa']) + "   ")
          fh.write('\n')
      except:
        pass
      fh.write("P: " + '\n')
      for j in range(len(g.top_parameters[i][1])):
        fh.write('{:12.4e}'.format(g.top_parameters[i][1][j]) + " ")
        if(j % 5 == 4):
          fh.write('\n')
      fh.write('\n')
      fh.write('\n')
    fh.close()


    potential.save(pot_dir, p)

    pf_save.save_summary(pot_dir, p, rss, bp, bp_known, efs, efs_known)

    b_props.bp_eos_plot(plot_dir)




  def save_summary(dir_out, p, rss, bp, bp_known, efs, efs_known):

    fh = open(dir_out + 'summary.txt', 'w')
    fh.write("RSS:  " + '\n')
    fh.write(str(rss) + '\n')
    fh.write("" + '\n')
    fh.write("P:  " + '\n')
    for pn in range(len(p)):
      fh.write(str(p[pn]) + '\n')
    fh.write("" + '\n')

    for i in range(len(bp_known['bp_calculations'])):  
      b = bp_known['bp_calculations'][i]   
      fh.write("Known:" + '\n')

      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(b['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(b['a0']) + '\n')
      fh.write("e0:     " + str(b['e0']) + '\n')
      fh.write("B0:     " + str(b['b0_gpa']) + '\n')

      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write('\n')


    for i in range(len(bp['bp_calculations'])):  
      b = bp['bp_calculations'][i]   
      fh.write("Calculated:" + '\n')

      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(b['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(b['a0']) + '\n')
      fh.write("e0:     " + str(b['e0']) + '\n')
      fh.write("B0:     " + str(b['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write("" + '\n')



    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.close()   


  def rss_plot(name, d):     
    d = numpy.asarray(d)
    dir_out = g.dirs['wd'] + '/fitting/saved_potentials/' + pf_save.count_name + '/'
    plot_dir = dir_out + 'rss_plot/'
    
    std.make_dir(dir_out)
    std.make_dir(plot_dir)

    plt.figure(figsize=(8,6))
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')   
    plt.title("RSS vs Time " + name)
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    plt.plot(d[:,0], d[:,1], 'k')
    plt.savefig(plot_dir + 'rss_over_time.eps', type='efs')
    plt.close('all') 

    plt.figure(figsize=(8,6))
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')   
    plt.title("RSS vs Time " + name)
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    plt.ylim(1.0,10**(numpy.log10(1.0*max(d[:,1])) + 1))
    plt.plot(d[:,0], d[:,1], 'k')
    plt.yscale("log")
    plt.savefig(plot_dir + 'rss_over_time_log.eps', type='efs')
    plt.close('all') 



  def rss_plot_full():     

    plot_dir = g.dirs['wd'] + '/fitting/rss_plot/'
    std.make_dir(plot_dir)

    plt.figure(figsize=(12,8))
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')  
    plt.title("RSS vs Time")
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    c = ['b','g','r','k']
    n = 0
    for p in pf.rss_plot_data:
      dt = p[0] - pf.start_time
      rss_name = p[1]
      xy = numpy.copy(p[2])
      xy[:,0] = xy[:,0] + dt
      plt.plot(xy[:,0], xy[:,1], c[n%len(c)], label=rss_name)
      n = n + 1
    plt.legend()
    plt.savefig(plot_dir + 'rss_over_time.eps', type='efs')
    plt.close('all') 


    plt.figure(figsize=(12,8))
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')  
    plt.title("RSS vs Time")
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    c = ['b','g','r','k']
    n = 0
    for p in pf.rss_plot_data:
      dt = p[0] - pf.start_time
      rss_name = p[1]
      xy = numpy.copy(p[2])
      xy[:,0] = xy[:,0] + dt
      plt.plot(xy[:,0], xy[:,1], c[n%len(c)], label=rss_name)
      n = n + 1
    plt.legend()
    plt.yscale("log")
    plt.savefig(plot_dir + 'rss_over_time_log.eps', type='efs')
    plt.close('all') 


    """
 
    plt.title("RSS vs Time " + name)
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    plt.plot(d[:,0], d[:,1], 'k')
    plt.savefig(plot_dir + 'rss_over_time.eps', type='efs')
    plt.close('all') 

    plt.figure(figsize=(8,6))
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')   
    plt.title("RSS vs Time " + name)
    plt.xlabel("Time/s")
    plt.ylabel("RSS")
    plt.ylim(1.0,10**(numpy.log10(1.0*max(d[:,1])) + 1))
    plt.plot(d[:,0], d[:,1], 'k')
    plt.yscale("log")
    plt.savefig(plot_dir + 'rss_over_time_log.eps', type='efs')
    plt.close('all') 
    """


   


    """  

    fh = open(dir_out + 'summary.txt', 'w')
    fh.write("RSS:  " + '\n')
    fh.write(str(g.pfdata['rss']['best']) + '\n')
    fh.write("" + '\n')
    fh.write("P:  " + '\n')
    for pn in range(len(g.top_parameters[0][1])):
      fh.write(str(g.top_parameters[0][1][pn]) + '\n')
    fh.write("" + '\n')


    for i in range(len(g.top_parameters[0][3]['bp_calculations'])):
      
      fh.write("Known:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.bp_known['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.bp_known['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.bp_known['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.bp_known['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write('\n')

      fh.write("Calculated:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.top_parameters[0][3]['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write("" + '\n')


    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.write("" + '\n')


    fh.write("RSS Details:  " + '\n')
    fh.write("===============================  " + '\n')
    fh.write("Total:      " + str(g.top_parameters[0][0]) + '\n')
    fh.write("a0:         " + str(bp.rss_by_type[0]) + '\n')
    fh.write("e0:         " + str(bp.rss_by_type[1]) + '\n')
    fh.write("b0:         " + str(bp.rss_by_type[2]) + '\n')
    fh.write("ec:         " + str(bp.rss_by_type[3]) + '\n')
    fh.write("g:          " + str(bp.rss_by_type[4]) + '\n')
    fh.write("e:          " + str(bp.rss_by_type[5]) + '\n')
    fh.write("v:          " + str(bp.rss_by_type[6]) + '\n')
    fh.write("a0 (w):     " + str(bp.rss_by_type_w[0]) + '\n')
    fh.write("e0 (w):     " + str(bp.rss_by_type_w[1]) + '\n')
    fh.write("b0 (w):     " + str(bp.rss_by_type_w[2]) + '\n')
    fh.write("ec (w):     " + str(bp.rss_by_type_w[3]) + '\n')
    fh.write("g (w):      " + str(bp.rss_by_type_w[4]) + '\n')
    fh.write("e (w):      " + str(bp.rss_by_type_w[5]) + '\n')
    fh.write("v (w):      " + str(bp.rss_by_type_w[6]) + '\n')


    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.close()   

    #potential.output(pot_dir, g.pfdata['p']['best'])
    """
    """
    p = numpy.copy(g.top_parameters[0][1])
    potential.update(p)
    rss = pf.get_rss()

    fh = open(dir_out + 'summary.txt', 'w')
    fh.write("RSS:  " + '\n')
    fh.write(str(g.top_parameters[0][0]) + '\n')
    fh.write("" + '\n')
    fh.write("P:  " + '\n')
    for pn in p:
      fh.write(str(pn) + '\n')
    fh.write("" + '\n')
    fh.write("EFS Density:" + '\n')
    fh.write(str(g.top_parameters[0][4]) + '\n')
    fh.write("" + '\n')
    fh.write("BP Density:" + '\n')
    fh.write(str(g.top_parameters[0][5]) + '\n')
    fh.write("" + '\n')
    fh.write("" + '\n')
    for i in range(len(g.top_parameters[0][3]['bp_calculations'])):
      
      fh.write("Known:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.bp_known['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.bp_known['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.bp_known['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.bp_known['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write('\n')

      fh.write("Calculated:" + '\n')
      ec = []
      for ei in range(6):
        ec.append([])
        for ej in range(6):
          ec[ei].append(str('{:6.3f}'.format(g.top_parameters[0][3]['bp_calculations'][i]['ec_gpa'][ei, ej])))

      fh.write("a0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['a0']) + '\n')
      fh.write("e0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['e0']) + '\n')
      fh.write("B0:     " + str(g.top_parameters[0][3]['bp_calculations'][i]['b0_gpa']) + '\n')
      for ei in range(6):
        if(ei == 0):          
          fh.write("EC:     ")
        else:          
          fh.write("        ")
        for ej in range(6):
          fh.write(ec[ei][ej] + ' ')
        fh.write('\n')
      fh.write("" + '\n')


    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.write("" + '\n')


    fh.write("RSS Details:  " + '\n')
    fh.write("===============================  " + '\n')
    fh.write("Total:      " + str(g.top_parameters[0][0]) + '\n')
    fh.write("a0:         " + str(bp.rss_by_type[0]) + '\n')
    fh.write("e0:         " + str(bp.rss_by_type[1]) + '\n')
    fh.write("b0:         " + str(bp.rss_by_type[2]) + '\n')
    fh.write("ec:         " + str(bp.rss_by_type[3]) + '\n')
    fh.write("g:          " + str(bp.rss_by_type[4]) + '\n')
    fh.write("e:          " + str(bp.rss_by_type[5]) + '\n')
    fh.write("v:          " + str(bp.rss_by_type[6]) + '\n')
    fh.write("a0 (w):     " + str(bp.rss_by_type_w[0]) + '\n')
    fh.write("e0 (w):     " + str(bp.rss_by_type_w[1]) + '\n')
    fh.write("b0 (w):     " + str(bp.rss_by_type_w[2]) + '\n')
    fh.write("ec (w):     " + str(bp.rss_by_type_w[3]) + '\n')
    fh.write("g (w):      " + str(bp.rss_by_type_w[4]) + '\n')
    fh.write("e (w):      " + str(bp.rss_by_type_w[5]) + '\n')
    fh.write("v (w):      " + str(bp.rss_by_type_w[6]) + '\n')


    fh.write("" + '\n')
    fh.write("" + '\n')
    fh.close()

    #potential.plot_python_potentials(plot_dir)
    potential.output(pot_dir)

    pf_top.save(dir_out, 'top_parameters.txt')
    """





