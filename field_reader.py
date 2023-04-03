# based on code made by collaborator Luke Bhan

from os import listdir
from os.path import isfile, join
import math

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
mypath = "/shared/home/malavjm1/hybrid_system/s0/better_gs/run_data"
prefix = "Eind_avgY_"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and prefix in f]
onlyfiles.sort()
n=0
writer = open('neg_field', 'w')
for idx, f in enumerate(onlyfiles):
    print(f)
    with open(f) as fp:
        line = fp.readline()
        sum = 0
        while line:
            arr = line.split()
            if(float(arr[0])) <= 3.557 and float(arr[0]) >= -11.20 and float(arr[1]) <=-16.02 and float(arr[1]) >=-48.30: #and float(arr[0]) > 0.0:
                print("TEST1")
                sum += float(arr[2])
                n=n+1

            line = fp.readline()

        sum = sum/n
        writer.write(str(float(remove_prefix(f, prefix))*0.02) + " " + str((sum))+ "\n")

print(onlyfiles)
#X_MIN_GRID=     -32.2009460
#X_MAX_GRID=      31.7474115
#Y_MIN_GRID=     -18.5949125
#Y_MAX_GRID=      18.1413780
#Z_MIN_GRID=     -48.3014190
#Z_MAX_GRID=      47.8478845
