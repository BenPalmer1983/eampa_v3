######################################################
#  Ben Palmer University of Birmingham 2020
#  Free to use
######################################################




class pgrad:

  def run():  
  
    print("Parameter Gradient") 

        
    # Setup EFS
    efs.init()                         # Initialise (allocate arrays)
    efs_calc.set_weights()
    potential.load_to_efs()
    configs.efs_add_config()           # Add configs

    # Setup BP
    bp_calc.init()
    bp_calc.set_weights()
    potential.load_to_bp()
    b_props.bp_add()


    # Run EFS and BP
    rss = rss_calc.run_calc()

    p = numpy.copy(potential.pot['p'])


    a = 0
    for fn in range(len(potential.pot['functions'])):
      print()
      print(potential.pot['functions'][fn]['file'])
      print("==============================================")
      print('{:3s}'.format(""), end=" ")
      print('{:17s}'.format("     RSS B"), end=" ")
      print('{:17s}'.format("     RSS"), end=" ")
      print('{:17s}'.format("     RSS F"), end=" ")
      print('{:17s}'.format("     PARAM B"), end=" ")
      print('{:17s}'.format("     PARAM"), end=" ")
      print('{:17s}'.format("     PARAM F"), end=" ")
      print('{:17s}'.format("     GRAD"), end=" ")
      print('{:17s}'.format("     PFIXED"), end=" ")
      print()
      if(potential.pot['functions'][fn]['f_on'] == 1):
        b = a + len(potential.pot['functions'][fn]['params'])
        c = 0
        for i in range(a, b,1):
          if(p[i] == 0.0):
            dh =1.0E-8
          else:
            dh = min(1.0E-8,1.0E-6 * p[i])

          # Forward
          p_f = numpy.copy(p)
          p_f[i] = p[i] + dh
          potential.update(p_f)
          rss_f = rss_calc.run_calc()

          # Backward
          p_b = numpy.copy(p)
          p_b[i] = p[i] - dh
          potential.update(p_b)
          rss_b = rss_calc.run_calc()


          p_grad = (rss_f - rss_b) / (2.0 * dh)

          print('{:3n}'.format(i), end=" ")
          print('{:17.9e}'.format(rss_b), end=" ")
          print('{:17.9e}'.format(rss), end=" ")
          print('{:17.9e}'.format(rss_f), end=" ")
          print('{:17.9e}'.format(p_b[i]), end=" ")
          print('{:17.9e}'.format(p[i]), end=" ")
          print('{:17.9e}'.format(p_f[i]), end=" ")
          print('{:17.9e}'.format(p_grad), end=" ")

          if(len(potential.pot['functions'][fn]['params_fixed'])>= len(potential.pot['functions'][fn]['params'])):
            print('{:17.9e}'.format(potential.pot['functions'][fn]['params_fixed'][c]), end=" ")
          print()
          c = c + 1
        a = b

    
    potential.update(p)





