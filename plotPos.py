import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

posX = []
posY = []
posZ = []

f = open('starting_configuration5.dat','r')
x, y, z = np.loadtxt(f, delimiter=' ', usecols=(0,1,2), unpack=True)

g = open('final_configuration5.dat','r')
xf, yf, zf = np.loadtxt(g, delimiter=' ', usecols=(0,1,2), unpack=True)

h = open('intermediate_configuration5.dat','r')
xm, ym, zm = np.loadtxt(h, delimiter=' ', usecols=(0,1,2), unpack=True)

i = open('intermediate_configuration2_5.dat','r')
xm2, ym2, zm2 = np.loadtxt(i, delimiter=' ', usecols=(0,1,2), unpack=True)

#Plot Vertices and Fist Hit
fig = plt.figure(1)
ax = Axes3D(fig)
#ax.plot(x,y,z,'ro', label='start')
#ax.plot(xm,ym,zm,'go', label='mid')
ax.plot(xm2,ym2,zm2,'co', label='mid2')
ax.plot(xf,yf,zf,'bo', label='final')
ax.legend(loc='upper left', numpoints = 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()