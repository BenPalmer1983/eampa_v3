

class pf_saved:

  def run():
  
    #pf_pgradient.run()
    g.pfdata['stage'] = 'START - LOAD SAVED'

    # Load saved parameters
    try:
      params = pf_top.load(g.dirs['wd'] + '/save', 'top_parameters.txt')    
      for p in params:
        pf_potential.update(p)
        rss = pf.get_rss()
    except:
      pass