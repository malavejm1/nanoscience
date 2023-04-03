# usual imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os

# for turning visual output into png
fg = plt.figure()

# names + locations of key files
fname="../grid_info.dat"
fname2="Eind_avgY_*"
fname3="../../t2/run_data/"+fname2

# initialize x-z grid points
input=np.genfromtxt(fname,dtype='str')
l=input[:,0]  #strings
n=input[:,1]  #strings to be used with float()


# get grid points in the x/z dirs
xgrid=[];zgrid=[];j=0
for ll in l:
 if ll == 'X':
  xgrid.append(float(n[j]))
 if ll == 'Z':
  zgrid.append(float(n[j]))
 j=j+1



# get list of file names
os.system("ls "+fname2+" > densfiles")
os.system("wc -l densfiles | \
 		    awk '{print $1}' > num")

# get number of files to be inspected
num_fname1=np.loadtxt("num",unpack=True)

# put names of all density filesi in array
file_bin1=np.genfromtxt("densfiles",dtype='str')
os.system("rm num densfiles")

os.system("ls "+fname3+" > densfiles")
os.system("wc -l densfiles | \
                    awk '{print $1}' > num")
num_fname2=np.loadtxt("num",unpack=True)
file_bin2=np.genfromtxt("densfiles",dtype='str')
os.system("rm num densfiles")

def plotter(iter):
# get and organize density data -----


# data from field file
 os.system("cat "+file_bin1[iter]+ \
		   " | awk '{print $5}' > data")
 m1 = np.loadtxt("data",unpack=True)
 os.system("rm data")

 os.system("cat "+file_bin2[iter]+ \
                   " | awk '{print $5}' > data")                                                                
 m2 = np.loadtxt("data",unpack=True)
 os.system("rm data")

# proper spatial organization of data
 data1 = np.reshape(m1,(len(xgrid),len(zgrid))
)
 data2 = np.reshape(m2,(len(xgrid),len(zgrid))
)
 


# generation of plots
 p1 = plt.subplot(211)
 c1 = plt.pcolormesh(zgrid,xgrid,data1**2,
			cmap='winter',
			shading='gouraud',
			vmin=1e-5,
			vmax=0.1)


 p2 = plt.subplot(212)
 c2 = plt.pcolormesh(zgrid,xgrid,data2**2,
                        cmap='winter',
                        shading='gouraud',
                        vmin=-0.1,
                        vmax=0.1)
 

# bells and whistles
 p1.set_title("HNP induced field: Al stair alignment")
 p1.set_xlabel("z")
 p1.set_ylabel("x")
 p2.set_title("HNP induced field: no Al stair alignment") 
 p2.set_xlabel("z")
 p2.set_ylabel("x")
 p1.annotate('', 
	xy=(-8,50), 
	xytext=(-10,50),
	arrowprops=dict(color='green',
				shrink=0.005,
				width=0.25,
				headwidth=5)
)


 plt.subplots_adjust(bottom=0.1,
	right=0.8,
	top=0.9,
	hspace=0.5)
 
 cax = plt.axes([0.85,0.1,0.075,0.8])
 plt.colorbar(cax=cax)
 
 fg.savefig("field"+str(iter)+".png")
 print("field"+str(iter)+".png done")
 plt.clf()
 return


# execution of user-defined function
if num_fname2 > num_fname1:
 num = int(num_fname1)
else:
 num = int(num_fname2)
for i in np.arange(0,num):
 plotter(i)
