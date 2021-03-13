
class pf_ng:

  R = None
  rss = None

  def run():
    g.pfdata['stage'] = 'NEWTON GAUSS'
    g.pfdata['stage_brief'] = 'NG'

    p = g.top_parameters[0][1]
    pf_ng.set(p)
    pf_ng.residual()

    Jw = len(p)
    Jh = len(pf_ng.R)
    R = numpy.copy(g.rss['residual'])
    h = 1.0e-10
    l = 0.1

    J = numpy.zeros((Jh,Jw),)
    for i in range(len(p)):
      p_b = numpy.copy(p)
      p_b[i] = p_b[i] - h
      pf_potential.update(p_b)      
      rss_b = pf.get_rss(False)
      R_b = numpy.copy(g.rss['residual'])

      p_f = numpy.copy(p)
      p_f[i] = p_f[i] + h
      pf_potential.update(p_f)      
      rss_f = pf.get_rss(False)
      R_f = numpy.copy(g.rss['residual'])

      J[:,i] = (R_f - R_b) / (2.0 * h)

    JT = numpy.transpose(J)
    JTJ = numpy.matmul(JT, J)
    JTJ_diag = numpy.diag(JTJ)
    JTR = numpy.matmul(JT, R) 

    try:
      dp = numpy.linalg.solve(JTJ + l * JTJ_diag, -JTR)
      p = p - dp
    except:
      pass
 

    pf_potential.update(p)
    pf.get_rss(False)



  def residual():
    rss = pf.get_rss(False)
    pf_ng.R = g.rss['residual']
    pf_ng.rss = rss
    

  def jacobian():
    print("J")







  def set(p):
    pf_potential.update(p)
