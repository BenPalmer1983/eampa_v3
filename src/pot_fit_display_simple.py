
class pf_display_simple:


  def output():
    if(g.pfdata['last_stage'] == None or g.pfdata['last_stage'] != g.pfdata['stage']):
      g.pfdata['last_stage'] = g.pfdata['stage']
      print("##########################################################################")
      print(g.pfdata['stage'])
      print("##########################################################################")
    
    rss_current = "Error"
    a0 = ''
    e0 = ''
    b0 = ''

    if(g.pfdata['rss']['current'] != None):
      rss_current = '{:12.4e}'.format(g.pfdata['rss']['current'])

    bplen = 14
    if(g.pfdata['bp']['current'] != None):
      bplen = 0 
      for i in range(len(g.pfdata['bp']['current']['bp_calculations'])):
        a0 = a0 + '{:12.4e}'.format(g.pfdata['bp']['current']['bp_calculations'][i]['a0']) + ' '
        e0 = e0 + '{:12.4e}'.format(g.pfdata['bp']['current']['bp_calculations'][i]['e0']) + ' '
        b0 = b0 + '{:12.4e}'.format(g.pfdata['bp']['current']['bp_calculations'][i]['b0_gpa']) + ' '
        bplen = bplen + 14
    

    pf_time = '{:6.3f}'.format(time.time() - pf.start_time)
    rss_best = '{:12.4e}'.format(g.pfdata['rss']['best'])

    sb = g.pfdata['stage_brief']
    while(len(sb)<8):
      sb = sb + " "

    print(sb, end = " ")
    print(pf_display.pad_l(pf_time, 10), end="    ")
    print(pf_display.pad_r(rss_best, 14), end=" ")
    print(pf_display.pad_r(rss_current, 14), end=" ")
    print(pf_display.pad_r("", 14), end=" ")
    print(pf_display.pad_r(a0, bplen), end=" ")
    print(pf_display.pad_r(e0, bplen), end=" ")
    print(pf_display.pad_r(b0, bplen), end=" ")
    print()




