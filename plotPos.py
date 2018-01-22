''' Plot Particle Configuration '''

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

posX = []
posY = []
posZ = []


# read files
name = 'P5_'

f = open(name + 'starting_configuration.dat','r')
x, y, z = np.loadtxt(f, delimiter=' ', usecols=(0,1,2), unpack=True)

g = open(name + 'final_configuration.dat','r')
xf, yf, zf = np.loadtxt(g, delimiter=' ', usecols=(0,1,2), unpack=True)

h = open(name + 'intermediate_configuration.dat','r')
xm, ym, zm = np.loadtxt(h, delimiter=' ', usecols=(0,1,2), unpack=True)

i = open(name + 'intermediate_configuration2.dat','r')
xm2, ym2, zm2 = np.loadtxt(i, delimiter=' ', usecols=(0,1,2), unpack=True)

j = open(name + 'intermediate_configuration3.dat','r')
xm3, ym3, zm3 = np.loadtxt(j, delimiter=' ', usecols=(0,1,2), unpack=True)

k = open(name + 'intermediate_configuration4.dat','r')
xm4, ym4, zm4 = np.loadtxt(k, delimiter=' ', usecols=(0,1,2), unpack=True)


#Plot Particle Positions
fig = plt.figure(1)
ax = Axes3D(fig)
#ax.plot(x,y,z,'ro', label='start')
#ax.plot(xm,ym,zm,'go', label='mid')
#ax.plot(xm2,ym2,zm2,'co', label='mid2')
#ax.plot(xm3,ym3,zm3,'mo', label='mid3')
#ax.plot(xm4,ym4,zm4,'ko', label='mid4')
ax.plot(xf,yf,zf,'bo', label='final')
ax.legend(loc='upper left', numpoints = 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()