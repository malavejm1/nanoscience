# a code that outputs a visual of where a super cell
# is located with respect to imaged cells in a system
# described using periodic boundary conditions 

import numpy as np
import matplotlib.pyplot as plt
import os

fg = plt.figure()


N_x=300
N_y=49
N_z=85
dx=0.474

ay = N_y*1.5
az = N_z*1.5


os.system('tail -n +2 atom.inp > animus')
atomx,atomy,atomz,s,m = np.loadtxt('animus',unpack=True)
os.system("rm animus")


xmax =-(N_x*0.5*dx)
xmin =(N_x*0.5*dx)
ymin=-(ay*1.5)
ymax=(ay*1.5)
zmin=-(az*1.5)
zmax=(az*1.5)


plt.plot(atomy,atomz,'o', color = 'red',alpha = 0.08)


miy = atomy - N_y/2
miz = atomz - N_z/2
piy = atomy + N_y/2
piz = atomz + N_z/2

plt.plot(miy,atomz,'o', color = 'green',alpha = 0.08)
plt.plot(piy,atomz,'o', color = 'green',alpha = 0.08)
plt.plot(atomy,miz,'o', color = 'green',alpha = 0.08)
plt.plot(atomy,piz,'o', color = 'green',alpha = 0.08)


plt.xlim(-50,50)
plt.ylim(-100,100)
fg.savefig('pbc.png')

