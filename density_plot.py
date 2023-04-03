#Strings together N data files to make N density plots or contour plots. I usually follow this up with ffmpeg, a way to string the output together to make a movie

# DENSITY PLOT PYTHON PROGRAM
# READ MY COMMENTS CAREFULLY


# the usual imports. you can ignore these
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib as mpl
import numpy as np
import pylab as py
import math as m
import os
fg = plt.figure()




x= 150   # grid points (put them here)
y = 16
z = 16




# user-defined function that sets up the plotting system. Try not to edit more
# than is necessary here

def plotter(iter,picture):
 ID = "{:07}".format(iter)
 pavx,pavz,cx,cy,cavz = np.loadtxt('curr_avgY_'+ID,unpack=True)
 px,pz,cz = np.loadtxt('dens_sliceY0_'+ID,unpack=True) # read in density slice file with certain ID number
 

 os.system("tail -" + str(x) + " 1d_Ezinc_Egrid" +ID + "  > file") # reads in laser file (ignore)
 os.system("awk '{print $2}' file > file2")

 print(ID)
 E0 = np.loadtxt('file2',unpack=True) # ignore
 os.system("rm file file2")

 dummy1,dummy2,density = np.loadtxt('dens_sliceY0_0000100',unpack=True)  # ignore this
 os.system("awk '{print $2}' grid_info.dat > GRID")
 grid = np.loadtxt("GRID",unpack=True); os.system("rm GRID")
 gridx = np.zeros(len(grid[0:x]),float);gridx = grid[0:x]              # ignore all of this
 gridy = np.zeros(len(grid[x:x+y]),float);gridy = grid[x:x+y]
 gridz = np.zeros(len(grid[x+y:x+y+z]),float);gridz = grid[x+y:x+y+z]

 #----------------------------------------------------------------------------------------------
 # READ IN ATOMIC COORDINATES
 # atomx, atomy, atomz, ditch1, ditch2 ---> all arrays that store the values in atom.inp
  # you have to make sure the atom.inp file is commensurate with this
 os.system('tail -n +2 "atom.inp" > animus')
 atomx,atomy,atomz, ditch1, ditch2 = np.loadtxt('animus',unpack=True)  # you may need to edit this
 os.system("rm animus")


 # RESOLUTION DIAL AND SAMPLING ALGORITHIM
 # for the love of god, be careful from here on.

 nf= np.zeros(len(px),float)
 for i in range(0,len(px)):
  #nf[i] = np.sqrt(cx[i]**2 + cy[i]**2 + cz[i]**2) #vector magnitude
  nf[i] = cz[i]

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
 Res = 7  # ignore


 atx,atz = np.zeros(len(px[slice(1,x,Res)])),np.zeros(len(px[slice(1,x,Res)]))
 ctx,ctz = np.zeros(len(px[slice(1,x,Res)])),np.zeros(len(px[slice(1,x,Res)]))
 tx = []
 tz = []
 tcx = []
 tcz = []

# SAMPLING ALGORITHM

# ignore

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


 #qmn = str(c[tab])



 c = plt.pcolormesh(gridx,gridz,field_matrix,cmap = "bwr") # you may need to edit this.
 # but you cannot usefully edit this until you plot the files
 # data from those files can be put here to refine the density plot visuals
 plt.colorbar(c)


 plt.plot(atomx,atomz,'o', color = 'black',alpha = 0.08)

 print("tab=",tab)
# if picture == 1: 
# py.ylabel("z (a.u.)")
# py.xlabel("x (a.u.)")
 py.title("Current density of graphene sheet (kick along x)")  # name of graph. change this
 py.ylabel("z (a.u.)")  #keep
 py.xlabel("x (a.u.)")  # keep
 fg.savefig('image'+str(tab)+'.png')

 print(ID + " finished")
 
 plt.clf()


 return


# images for movie made here
# this part makes the images

stride = 100   #put the three digit number for the density file here. I think it was '312' for you
journey = 5000 # just make this some big number. Bigger than the number of your final density file
tab = 1
c = np.zeros(int(journey/stride),int)
for k in range(stride,journey+stride,stride):
 plotter(k,1)
 tab = tab+1
