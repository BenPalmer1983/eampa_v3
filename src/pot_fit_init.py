class pf_init:

  def run():  

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
    bp_calc.get_known()

    g.pfdata = {}
    
    g.pfdata['stage'] = ''
    g.pfdata['last_stage'] = None

    # Set up RSS
    g.pfdata['rss'] = {'start': None,
                       'current': None,
                       'best': None,
                       'counter': 0,
                       'counter_successful': 0,
                       'since_improvement': 0,}


    g.pfdata['psize'] = potential.pot['p_count']
    g.pfdata['p'] = {'current': None, 'best': None, 'known': None,}
    g.pfdata['bp'] = {'current': None, 'best': None, 'known': None,}
    g.pfdata['efs'] = {'current': None, 'best': None,}
    #g.pfdata['max_density'] = {'bp_current': None, 'efs_current': None, 'bp_best': None, 'efs_best': None,}


    # SAVE Starting Parameters & Variation
    g.pfdata['params_start']  = numpy.copy(potential.pot['p'])
    g.pfdata['params_var'] = numpy.copy(potential.pot['p_var'])



    g.pfdata['bp']['known'] = g.bp_known.copy()
    g.pfdata['efs']['known'] = g.efs_known.copy()












