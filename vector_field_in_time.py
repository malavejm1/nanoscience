
# vector plot of current in nanostructure
# Shows the flow of the electrons/field in
# real-time. Very similar to density_plot.py


# the usual imports


import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib as mpl
import numpy as np
import pylab as py
import math as m
import os



fg = plt.figure()
x= 142   # grid points
y = 82
z = 213


px0,pz0,cx0,cy0,cz0 = np.loadtxt('curr_avgY_0000000',unpack=True)



def plotter(iter,picture):
 ID = "{:07}".format(iter)
 px,pz,cx,cy,cz = np.loadtxt('curr_avgY_'+ID,unpack=True)
 
# os.system("tail -" + str(x) + " 1d_Ezinc_Egrid" +ID + "  > file")
# os.system("awk '{print $2}' file > file2")

 print(ID)
# E0 = np.loadtxt('file2',unpack=True)
# os.system("rm file file2")

 dummy1,dummy2,density = np.loadtxt('GSOD_avgY_0000000',unpack=True)
 os.system("awk '{print $2}' ../grid_info.dat > GRID")
 grid = np.loadtxt("GRID",unpack=True); os.system("rm GRID")
 gridx = np.zeros(len(grid[0:x]),float);gridx = grid[0:x] 
 gridy = np.zeros(len(grid[x:x+y]),float);gridy = grid[x:x+y]
 gridz = np.zeros(len(grid[x+y:x+y+z]),float);gridz = grid[x+y:x+y+z]

 #----------------------------------------------------------------------------------------------
 # READ IN ATOMIC COORDINATES
 os.system('tail -n +2 ../atom.inp > animus')
 atomx,atomy,atomz,vvv,bbb = np.loadtxt('animus',unpack=True)
 os.system("rm animus")
 # RESOLUTION DIAL AND SAMPLING ALGORITHIM


 nf= np.zeros(len(px),float)
 for i in range(0,len(px)):
  nf[i] = cz[i]-cz0[i] #vector magnitude

 field_matrix = np.zeros([z,x],float)
 add = 0

# cutoff = 0.5*((np.max(nf) + np.min(nf)) / (np.max(E0) + np.min(E0)))

 for i in range(0,x):
  for j in range(0,z):
   field_matrix[j][i] = nf[j + add]
  add = add + z


 print(field_matrix.size)
 print(field_matrix.shape)




# must be an integer
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
   tcz.append(cz_entry)  #these arrays are used for vector plots


 #qmn = str(c[tab])



 plt.plot(atomx,atomz,'o',color = 'black',alpha = 0.08)
 X,Y = np.meshgrid(tx,tz)
 plt.quiver(X,Y,tcx,tcz)
 print(tab,X,Y)
# if picture == 1: 
# py.ylabel("z (a.u.)")
# py.xlabel("x (a.u.)")
 py.title("charge transfer: "+str(iter))
 py.ylabel("z (a.u.)")
 py.xlabel("x (a.u.)")  
 fg.savefig('image'+str(tab)+'.png')

 print(ID + " finished")
 
 plt.clf()


 return
# images for movie made here

stride = 100
journey = 25000
tab = 1
c = np.zeros(int(journey/stride),int)
for k in range(0,journey+stride,stride):
 plotter(k,1)
 tab = tab+1

