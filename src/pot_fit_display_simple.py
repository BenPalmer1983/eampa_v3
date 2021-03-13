
class pf_display_simple:


  def output():
    if(g.pfdata['last_stage'] == None or g.pfdata['last_stage'] != g.pfdata['stage']):
      g.pfdata['last_stage'] = g.pfdata['stage']
      print("#####################################")
      print(g.pfdata['stage'])
      print("#####################################")
    
    if(g.pfdata['rss']['current'] == None):
      rss_current = "Error"
    else:
      rss_current = '{:12.4e}'.format(g.pfdata['rss']['current'])

    pf_time = '{:6.3f}'.format(time.time() - pf.start_time)
    rss_best = '{:12.4e}'.format(g.pfdata['rss']['best'])

    sb = g.pfdata['stage_brief']
    while(len(sb)<8):
      sb = sb + " "

    print(sb, end = " ")
    print(pf_display.pad_l(pf_time, 10), end="    ")
    print(pf_display.pad_r(rss_best, 14), end=" ")
    print(pf_display.pad_r(rss_current, 14), end=" ")
    print()





