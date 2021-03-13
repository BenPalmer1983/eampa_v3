"""
Purely pseudo random 
Uses the upper and lower range as set in the potential fit files
Picks a random parameter centered on that range
Also uses the oversized_parameters setting to randomly pick an over sized range of parameters

e.g. fitting oversized_parameters=2.0,0.25,0.5,0.5,0.5

rn = numpy.random.uniform()    
b = g.fitting['oversized_parameters'][0]
prb = 1.0
  for i in range(1, len(g.fitting['oversized_parameters']),1):
    prb = prb * g.fitting['oversized_parameters'][i]
    if(rn<prb):
      m = m * b
    else:
      break

base is 2
if(rn >= 0.25) m = 1.0
if(rn < 0.25) m = 2.0
if(rn < 0.125) m = 4.0
if(rn < 0.0625) m = 8.0

"""

class pf_random:

  count = 0
  time_spent = 0.0

  def run():
    t_start = time.time()
    start_best_rss = g.pfdata['rss']['best'] 
    pf_random.count = pf_random.count + 1
    g.pfdata['stage'] = 'RANDOM ' + str(pf_random.count)
    g.pfdata['stage_brief'] = 'R' + str(pf_random.count)

    # START PARAMETER
    # Try randomly generated parameters
    for n in range(pf.fit['random_size']):
      p = pf_parameters.random_p(0.0)
      pf_potential.update(p)
      rss = pf.get_rss()

    pf_random.time_spent = pf_random.time_spent + (time.time() - t_start)

    pf.summary_line("RANDOM SEARCH " + str(pf_random.count), t_start, time.time(), start_best_rss,  g.pfdata['rss']['best'])
   
    pf_save.top("RANDOM_SEARCH_" + str(pf_random.count))


