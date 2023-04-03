import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


v1 = np.loadtxt('v1.dat')
v2 = np.loadtxt('v2.dat')

n = 101
X = np.zeros((n))
Y = np.zeros((n))
Z = np.zeros((n,n))
X = np.loadtxt('x.dat')
Y = np.loadtxt('y.dat')
Z = np.loadtxt('z.dat')
X,Y = np.meshgrid(X,Y)
Z=Z.reshape(n,n)


fig=plt.figure()
ax1=fig.add_subplot(111)

#fig, (ax0, ax1) = plt.subplots(2, 1)

c = ax1.pcolormesh(X, Y, Z, cmap=plt.cm.jet,vmin =Z.min(),vmax =Z.max(),shading='gouraud')
fig.colorbar(c, ax=ax1)

plt.xlim(-5,5)
plt.ylim(-5,5)
plt.show()
