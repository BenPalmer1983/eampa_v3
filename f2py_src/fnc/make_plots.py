
#export PYTHONPATH=$PYTHONPATH:"/cloud/Code/python/eampa/f2py_lib"


from f_fnc import fnc
from eampa_lib.f_spline import spline
import numpy
import matplotlib.pyplot as plt




class make_plots:


  
  @staticmethod
  def run():
  
    #####################################
    # PAIR
    #####################################
    
  
    x = numpy.linspace(1.0, 7, 100)    
    
    # LJ
    p = [2.3,3.5]
    p_fixed = [0.0]           # No fixed parameters
    y = fnc.lennard_jones_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('lennard_jones', 'Lennard-Jones', x, y) 
    
    # MORSE
    p = [4.669,1.256,2.8]
    p_fixed = [0.0]           # No fixed parameters
    y = fnc.morse_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('morse', 'Morse', x, y) 
    
    # BUCKINGHAM
    p = [6.0,0.5,12.0]
    p_fixed = [0.0]           # No fixed parameters
    y = fnc.buckingham_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('buckingham', 'Buckingham', x, y) 
    
    # QUARTIC POLYNOMIAL WITH REPULSIVE TERM
    p = [-4.5377,4.0659,-0.8548, 9.5272e7, 0.1193]
    p_fixed = [3.733]
    y = fnc.quartic_poly_rep_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('quartic_poly_rep', 'quartic_poly_rep', x, y) 
    
    
    
    
    #####################################
    # DENSITY
    #####################################
    
    x = numpy.linspace(0.0, 7, 100) 
    
    p = [3.816]
    p_fixed = [0.0]
    y = fnc.quadratic_density_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('quadratic_density', 'Quadratic Density', x, y)  
    
    
    p = [5.0, 1.323]
    p_fixed = [0.0]
    y = fnc.slater_4s_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('slater_4s', 'Slater 4s', x, y)  
    
    
    
    
    
  
    #####################################
    # EMBEDDING
    #####################################
            
    x = numpy.linspace(0.0, 1.0, 100) 
    
    p = [10.0]
    p_fixed = [0.0]
    y = fnc.fs_embedding_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('fs_embedding', 'Finnis-Sinclair Embedding', x, y)    
    
    p = [10.0]
    p_fixed = [0.0]
    y = fnc.mendelev_embedding_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('mendelev_embedding', 'Mendelev Embedding', x, y)      
    
    p = [-4., 3.5, -0.1]
    p_fixed = [0.0]
    y = fnc.triple_embedding_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('triple_embedding', 'Triple Embedding', x, y)       
    
    p = [-4.5, 2.95, -0.2]
    p_fixed = [0.0]
    y = fnc.ackland_embedding_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('ackland_embedding', 'Ackland Embedding', x, y)  
    
    
  
  
  
  
    #####################################
    # GENERAL
    #####################################
    
    x = numpy.linspace(0.0, 7.0, 100) 
    p = [0.9, 2.0, -0.1, 0.1, 0.07]
    p_fixed = [1.0, 2.0, 3.0, 5.0, 7.0]
    y = fnc.cubic_spline_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('cubic_spline', 'Cubic Spline', x, y)     
    
    x = numpy.linspace(0.0, 7.0, 100) 
    p = [0.9, 2.0, -0.1, 0.1, 0.07]
    p_fixed = [1.0, 2.0, 3.0, 5.0, 7.0]
    y = fnc.quintic_spline_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('quintic_spline', 'Quintic Spline', x, y)     
    
    
    nodes = numpy.asarray([[0.0,2.1],[1.0,0.5],[2.0,1.0],[3.0,0.2],[4.0,-0.5],[5.0,0.05],[7.0,0.0]])

    #xy = spline.get_nodes(x, y, 0.0, 7.0, 100, )
    xy = spline.spline_nodes(1, nodes, 100)
    make_plots.plot('poly3_node_spline', 'Poly3 Node Spline', xy[:,0], xy[:,1])   

    xy = spline.spline_nodes(2, nodes, 100)
    make_plots.plot('poly5_node_spline', 'Poly5 Node Spline', xy[:,0], xy[:,1])   


#S
#TYPE poly3
#X   0.0        0.65       1.3        1.95       2.6        3.25       3.9        4.55       5.2        5.85       6.5        
#Y   1.0        1.02352977 1.18823904 -0.0162325 0.04562567 0.02066891 0.00685261 0.00553381 0.00331582 0.00160056 0.0  
#POINTS 1001

    #x = numpy.linspace(0.0, 7.0, 100) 
    #p = [1000, 20, 1.0,-0.1,0.1,0.0]
    #p_fixed = [0.0, 1.0, 2.0, 3.0, 4.0, 7.0]
    #y = fnc.quintic_spline_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    #make_plots.plot('quintic_spline', 'Quintic Spline', x, y)  
    
    

    
    #spline_n_node(r, p, pf)
    
    
    
     
    
    """
    
    # Summed spline (Olsson et al)
    p = [0.976,-165.0,1.15,-78.499908, 1.216,-78.15495,1.650,1.8679553]
    p_fixed = [0.976,-165.0,1.15,-78.499908, 1.216,-78.15495,1.650,1.8679553]
    y = fnc.summed_spline_v(x, numpy.asarray(p, dtype=numpy.float32))
    make_plots.plot('summed_spline_pair', 'Summed Spline', x, y) 
    
    
    p = [0.344338030161, -10.0585662251, 2.19163816627, 2.71982098393, 17.415590781, 0.0020894685283]
    y = fnc.summed_spline_v(x, numpy.asarray(p, dtype=numpy.float32))
    make_plots.plot('summed_spline_pair', 'Summed Spline', x, y) 
    
    
    
    #P 0.344338030161 -10.0585662251 2.19163816627 2.71982098393 17.415590781 0.0020894685283
        
    p = [-4.5377,4.0659,-0.8548]
    p_fixed = [3.733]
    x = numpy.linspace(0.0, 7, 100)
    y = fnc.quartic_poly_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('quartic_poly', 'Quartic Polynomial', x, y)     
        
    p = [-4.5377,4.0659,-0.8548, 9.5272e7, 0.1193]
    p_fixed = [3.733]
    x = numpy.linspace(0.0, 7, 100)
    y = fnc.quartic_poly_rep_v(x, numpy.asarray(p, dtype=numpy.float32), numpy.asarray(p_fixed, dtype=numpy.float32))
    make_plots.plot('quartic_poly_rep', 'Quartic Polynomial With Repulsion', x, y)    
    
  
    #####################################
    # DENSITY
    #####################################
    
    
    p = [0.963,-11.0828,1.284,0.013905,1.685,-0.447541]
    y = fnc.summed_spline_v(x, numpy.asarray(p, dtype=numpy.float32))
    make_plots.plot('summed_spline_density', 'Summed Spline - Density', x, y)  
    
    """
    
    
    
    
    
    
    
    

  @staticmethod
  def plot(plot_file, plot_name, x, y, y_a = None, y_b = None):
    plt.clf()   
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(7,5))
    fig.tight_layout(pad=5.0)    
    fig.suptitle(plot_name)
    axs.plot(x, y, color='k', ls='solid', label='potential')
    if(type(y_a) == numpy.ndarray):
      axs.plot(x, y_a, color='k', ls='dashed', label='repulsion')
    if(type(y_b) == numpy.ndarray):
      axs.plot(x, y_b, color='k', ls='dotted', label='attraction')
    #axs.set_ylim(-5.0,5.0)
    axs.legend()
    plt.savefig('plots/' + plot_file + '.eps', format='eps')




make_plots.run()
