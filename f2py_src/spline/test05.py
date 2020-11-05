import numpy
from f_spline import spline
import matplotlib.pyplot as plt


nodes = numpy.zeros((6,2,),)
nodes[0,0] = 0.0
nodes[0,1] = 100.0
nodes[1,0] = 1.0
nodes[1,1] = 20.0
nodes[2,0] = 2.0
nodes[2,1] = -0.2
nodes[3,0] = 3.0
nodes[3,1] = 0.2
nodes[4,0] = 4.0
nodes[4,1] = -0.1
nodes[5,0] = 7.0
nodes[5,1] = 0.0

#([[0,1,2,3,7],[100,20,-0.2,0.2,0]])
print(nodes)

splined = spline.spline_nodes(2, nodes, 100)


print(splined)


plt.plot(nodes[:,0], nodes[:,1])
plt.plot(splined[:,0], splined[:,1])
plt.show()


