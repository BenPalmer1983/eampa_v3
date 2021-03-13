






class pf_save:
 
  count = 0

  def top(name):
    pf_save.count = pf_save.count + 1
 
    c = str(pf_save.count)
    while(len(c)<4):
      c = '0' + c
    name = c + '_' + name

    dir_out = g.dirs['wd'] + '/fitting/saved_potentials/' + name + '/'
    plot_dir = dir_out + 'plots/'
    pot_dir = dir_out + 'pots/'
    std.make_dir(dir_out)
    std.make_dir(plot_dir)
    std.make_dir(pot_dir)
    
    p = numpy.copy(g.top_parameters[0][1])
    pf_potential.update(p)
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

    potential.plot_python_potentials(plot_dir)
    potential_output.full(pot_dir)

    pf_top.save(dir_out, 'top_parameters.txt')





