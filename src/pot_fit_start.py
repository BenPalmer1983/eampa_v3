

class pf_start:

  def run():
  
    pf_display.stage = 'START - LOAD SAVED'
    
    p = potential.get_start_parameters()
    pf_potential.update(p)
    rss = pf.get_rss()

    #pf_pgradient.run()

    # Load saved parameters
    params = pf_top.load(g.dirs['wd'] + '/save', 'top_parameters.txt')    
    for p in params:
      pf_potential.update(p)
      rss = pf.get_rss()













