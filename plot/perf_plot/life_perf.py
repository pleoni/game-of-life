import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_acc.out"
datafile2="life_host.out"
datafile3="life_mic.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="dev_perf_eurora.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(211)
bottom = fig.add_subplot(212)

#############  TOP 

ACC_C1000 = data1[where((data1[:,0]==1) & (data1[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
ACC_C0    = data1[where((data1[:,0]==1) & (data1[:,5]==0)    ),:][0] # mpi 1 - comp 0

HOST_C1000 = data2[where((data2[:,0]==1) & (data2[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
HOST_C0    = data2[where((data2[:,0]==1) & (data2[:,5]==0)    ),:][0] # mpi 1 - comp 0

MIC_C1000 = data3[where((data3[:,0]==1) & (data3[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
MIC_C0    = data3[where((data3[:,0]==1) & (data3[:,5]==0)    ),:][0] # mpi 1 - comp 0

top.set_title(str(today) + ' life_hpc2 on eurora - NCOMP=1000')
top.grid()
top.set_xlabel('Lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

top.plot(ACC_C1000[:,3],ACC_C1000[:,8],'-xr',HOST_C1000[:,3],HOST_C1000[:,8],'-xg',MIC_C1000[:,3],MIC_C1000[:,8],'-xc');
top.legend(('KEP20','HOST','MIC'), loc = 'upper left', shadow = False, prop={'size':9})


#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc2 on eurora - NCOMP=0')
bottom.grid()
bottom.set_xlabel('Lattice size')
bottom.set_ylabel('time')

bottom.plot(ACC_C0[:,3],ACC_C0[:,8],'-xr',HOST_C0[:,3],HOST_C0[:,8],'-xg',MIC_C0[:,3],MIC_C0[:,8],'-xc');
bottom.legend(('KEP20','HOST','MIC'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


