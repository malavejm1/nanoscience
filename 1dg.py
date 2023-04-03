#code display spin-sensitive electron density in 1D
# along two, mutually orthogonal atomic/molecular axes

import numpy as np
import matplotlib.pyplot as m

tag = '5es2.5'
f1 = '/w0.5'
f2 ='/w1.0'
f3 ='/w1.5'

x1,dx1 = np.loadtxt(tag+f1+'/fort.11',unpack=True)
x2,dx2 = np.loadtxt(tag+f2+'/fort.11',unpack=True)
x3,dx3 = np.loadtxt(tag+f3+'/fort.11',unpack=True)
y1,dy1 = np.loadtxt(tag+f1+'/fort.21',unpack=True)
y2,dy2 = np.loadtxt(tag+f2+'/fort.21',unpack=True)
y3,dy3 = np.loadtxt(tag+f3+'/fort.21',unpack=True)

fig, axs = m.subplots(1,2,sharey=True)
axs[0].plot(x1,dx1,linestyle = 'dotted',color='red',label = "$\omega_{0}$ = 0.5")
axs[0].plot(x1,dx2,linestyle = 'dashed',color='blue',label = "$\omega_{0}$ = 1.0")
axs[0].plot(x1,dx3,linestyle = 'solid',color='purple',label = "$\omega_{0}$ = 1.5")
axs[0].set_xlim(-5,5)
axs[0].set_xlabel("x (a.u.)")
axs[0].set_ylabel("Spin-Up Density: $|\Phi|^2$")
axs[0].set_ylim(0,dx3.max()+.01)



axs[1].plot(y1,dy1,linestyle = 'dotted',color='red',label = "$\omega_{0}$ = 0.5")
axs[1].plot(y1,dy2,linestyle = 'dashed',color='blue',label = "$\omega_{0}$ = 1.0")
axs[1].plot(y1,dy3,linestyle = 'solid',color='purple',label = "$\omega_{0}$ = 1.5")
axs[1].set_xlim(-5,5)
axs[1].set_xlabel("y (a.u.)")
axs[1].set_ylim(0,dy3.max()+.01)

m.subplots_adjust(wspace=0,hspace=0)

m.legend(loc="upper center",ncol = 3,bbox_to_anchor=(0,1.1))

m.show()

