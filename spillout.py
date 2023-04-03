import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib as mpl
import numpy as np
import pylab as py
import math as m
import os
import matplotlib.colors as cc


def plotter(PATH,x,y,z):
 px,pz,cz = np.loadtxt(PATH+'/run_data/GSOD_sliceY0_0000001',unpack=True)
 cx = cz
 cy = cz

 #make grid
 os.system("awk '{print $2}' "+PATH+"/grid_info.dat > GRID")
 grid = np.loadtxt("GRID",unpack=True); os.system("rm GRID")
 gridx = np.zeros(len(grid[0:x]),float);gridx = grid[0:x]
 gridy = np.zeros(len(grid[x:x+y]),float);gridy = grid[x:x+y]
 gridz = np.zeros(len(grid[x+y:x+y+z]),float);gridz = grid[x+y:x+y+z]

# READ IN ATOMIC COORDINATES
 os.system('tail -n +2 '+PATH+'/atom.inp > animus')
 atomx,atomy,atomz,vvv,bbb = np.loadtxt('animus',unpack=True)
 os.system("rm animus")
 # RESOLUTION DIAL AND SAMPLING ALGORITHIM

 nf= np.zeros(len(px),float)
 for i in range(0,len(px)):
  nf[i] = cz[i] #vector magnitude

 field_matrix = np.zeros([z,x],float)
 add = 0

# cutoff = 0.5*((np.max(nf) + np.min(nf)) / (np.max(E0) + np.min(E0)))

 for i in range(0,x):
  for j in range(0,z):
   field_matrix[j][i] = nf[j + add]
  add = add + z


 Res = 7


 atx,atz = np.zeros(len(px[slice(1,x,Res)])),np.zeros(len(px[slice(1,x,Res)]))
 ctx,ctz = np.zeros(len(px[slice(1,x,Res)])),np.zeros(len(px[slice(1,x,Res)]))
 tx = []
 tz = []
 tcx = []
 tcz = []

# SAMPLING ALGORITHM

 for i in range(0,len(atx)):
  atx = px[slice(z + i*(Res - 1)*z - (z-1),z + i*(Res - 1)*z ,Res)]
  atz = pz[slice(z + i*(Res - 1)*z - (z-1),z + i*(Res - 1)*z ,Res)]
  ctx = cx[slice(z + i*(Res - 1)*z - (z-1),z + i*(Res - 1)*z ,Res)]
  ctz = cz[slice(z + i*(Res - 1)*z - (z-1),z + i*(Res - 1)*z ,Res)]
  for j in range(0,len(atx)):
   tx_entry,tz_entry = atx[j],atz[j]
   cx_entry,cz_entry = ctx[j],ctz[j]
   tx.append(tx_entry)
   tz.append(tz_entry)
   tcx.append(cx_entry)
   tcz.append(cz_entry)

 return atomx,atomz,gridx,gridz,field_matrix

x5nm,z5nm,gx5nm,gz5nm,f5nm = plotter('5nm',315,15,15)
x10nm,z10nm,gx10nm,gz10nm,f10nm = plotter('10nm',472,15,15)     
x20nm,z20nm,gx20nm,gz20nm,f20nm = plotter('20nm',864,15,15)     
x30nm,z30nm,gx30nm,gz30nm,f30nm = plotter('30nm',1218,15,15)     

tune1=(10**(15))*f5nm.min()
tune2=f5nm.max()/f30nm.max()


print(f5nm.min(),f5nm.max())
print(f10nm.min(),f10nm.max())
print(f20nm.min(),f20nm.max())
print(f30nm.min(),f30nm.max())


nm5min = f5nm.min()/tune1; nm5max = f5nm.max()/tune2
nm10min = f5nm.min()/tune1; nm10max = f5nm.max()/tune2
nm20min = f5nm.min()/tune1; nm20max = f5nm.max()/tune2
nm30min = f5nm.min()/tune1; nm30max = f5nm.max()/tune2

print(nm5min,nm5max)



fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=4)

cccc = ax1.pcolormesh(gx5nm,gz5nm,f5nm,norm=cc.LogNorm(vmin=nm5min, vmax = nm5max), cmap= "jet",shading='gouraud')
ax1.set_xlim(x5nm.min()-20,x5nm.min())
ax1.plot(x5nm,z5nm,'o', color = 'pink',alpha = 1)
#ax1.set_ylabel("z (a.u.)")
#ax1.set_xlabel("x (a.u.)")
ax1.text(x5nm.min()-19.5, 1, '5nm', fontsize=12,  color='white')

ax2.pcolormesh(gx10nm,gz10nm,f10nm,norm=cc.LogNorm(vmin=nm10min, vmax =nm10max ), cmap= "jet",shading='gouraud')
#plt.colorbar(nm5)
ax2.set_xlim(x10nm.min()-20,x10nm.min())
ax2.plot(x10nm,z10nm,'o', color = 'pink',alpha = 1)
#ax1.set_ylabel("z (a.u.)")
#ax1.set_xlabel("x (a.u.)")
ax2.text(x10nm.min()-19.5, 1, '10nm', fontsize=12,  color='white')


ax3.pcolormesh(gx20nm,gz20nm,f20nm,norm=cc.LogNorm(vmin=nm20min, vmax = nm20max), cmap= "jet",shading='gouraud')
#plt.colorbar(nm5)
ax3.set_xlim(x20nm.min()-20,x20nm.min())
ax3.plot(x20nm,z20nm,'o', color = 'pink',alpha = 1)
#ax1.set_ylabel("z (a.u.)")
#ax1.set_xlabel("x (a.u.)")
ax3.text(x20nm.min()-19.5, 1, '20nm', fontsize=12,  color='white')


spill = ['1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0.0']
tiks = np.linspace(x30nm.min()-20,x30nm.min(),11,endpoint=True)

v=ax4.pcolormesh(gx30nm,gz30nm,f30nm,norm=cc.LogNorm(vmin=nm30min, vmax = nm30max), cmap= "jet",shading='gouraud')
ax4.set_xlim(x30nm.min()-20,x30nm.min())
ax4.plot(x30nm,z30nm,'o', color = 'pink',alpha = 1)
#ax1.set_ylabel("z (a.u.)")
#ax1.set_xlabel("x (a.u.)")
ax4.set_xlabel("Spillout Distance (nm)")
ax4.text(x30nm.min()-19.5, 1, '30nm', fontsize=12,  color='white')




ax1.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
ax2.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
ax3.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
ax4.tick_params(left = False, right = False , labelleft = False)
ax4.set_xticks(tiks)
ax4.set_xticklabels(spill,minor=False)
plt.subplots_adjust(wspace=0, hspace=0)
fig.colorbar(cccc, ax=[ax1,ax2,ax3,ax4])
plt.show()
