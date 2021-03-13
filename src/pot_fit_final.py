
class pf_final:



  def run():


    p = g.pfdata['p']['best']
    pf_potential.update(p)

    # Run EFS and BP
    rss = rss_calc.run_calc()

   

    print('') 
    print('') 
    print('')     
    print('###########################################################################') 
    print('###########################################################################') 
    print('FINAL CALCULATION WITH BEST PARAMETERS') 
    print('###########################################################################') 
    print('###########################################################################') 
    print('')     


    b_props.bp_output_terminal()


    print('')     
    print('CONFIGS')    
    for n in range(efs.cc):    
      print('Config ' + str(n+1) + ':', efs.config_energy[n,2], efs.energies[n], (efs.config_energy[n,2]-efs.energies[n])**2)
    print('All configs:                   ' + str(efs.total_rss))
    print('All configs (energy):          ' + str(efs.energy_rss))
    print('All configs (force):           ' + str(efs.force_rss))
    print('All configs (stress):          ' + str(efs.stress_rss))
    print('All configs weighted:          ' + str(efs.total_rss_weighted))
    print('All configs weighted (energy): ' + str(efs.energy_rss_weighted))
    print('All configs weighted (force):  ' + str(efs.force_rss_weighted))
    print('All configs weighted (stress): ' + str(efs.stress_rss_weighted))

   
    print('')   
    for bp_id in range(bp.bp_configs_count):  
      print('')   
 
      print('BP  ' + str(bp_id)) 
      print('================') 

      rss_calc.print_line('a0', bp.known_alat[bp_id], bp.calc_alat[bp_id])
      rss_calc.print_line('v0', None, bp.calc_v0[bp_id])
      rss_calc.print_line('e0', bp.known_e0[bp_id], bp.calc_e0[bp_id])
      rss_calc.print_line('b0', bp.known_b0[bp_id], bp.calc_b0[bp_id])
      print('')   

      print("Calculated Stiffness Matrix (GPA)")
      for i in range(6):
        for j in range(6):
          print(str('{:14.6f}'.format(float(160.230732254e0 * bp.calc_ec[bp_id,i,j]))), end="")
        print()
      print("Known Stiffness Matrix (GPA)")
      for i in range(6):
        for j in range(6):
          print(str('{:14.6f}'.format(float(160.230732254e0 * bp.known_ec[bp_id,i,j]))), end="")
        print()
        #print(160.230732254e0 * bp.known_ec[bp_id,i,:])
      
      print('')   
  
    print('RSS:')    
    print('a0:', g.rss['bp']['a0'], "   w: ",g.rss['bp']['a0_weighted'])
    print('e0:', g.rss['bp']['e0'], "   w: ",g.rss['bp']['e0_weighted'])
    print('b0:', g.rss['bp']['b0'], "   w: ",g.rss['bp']['b0_weighted'])
    print('ec:', g.rss['bp']['ec'], "   w: ",g.rss['bp']['ec_weighted'])
    print('g:', g.rss['bp']['g'], "   w: ",g.rss['bp']['g_weighted'])
    print('e:', g.rss['bp']['e'], "   w: ",g.rss['bp']['e_weighted'])
    print('v:', g.rss['bp']['v'], "   w: ",g.rss['bp']['v_weighted'])
    print('')
    print('')
    print('RSS: ' + str(rss))
    print('')
    print('')
    print('Max Density:', bp.max_density, efs.max_density)
    print('')
    print('')
    #print(g.rss)
    print('')
    










    