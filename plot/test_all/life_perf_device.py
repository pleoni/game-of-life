import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_acc_perf.out"
datafile2="life_host_perf.out"
datafile3="life_mic_perf.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="device_perf_eurora.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(211)
bottom = fig.add_subplot(212)

#############  TOP 

ACC_C1000 = data1[where((data1[:,0]==1) & (data1[:,2]==16384) ),:][0] # mpi 1 - grid 16384
ACC_C0    = data1[where((data1[:,0]==1) & (data1[:,5]==1000)    ),:][0] # mpi 1 - comp 1000

HOST_C1000 = data2[where((data2[:,0]==1) & (data2[:,2]==16384) ),:][0] # mpi 1 - grid 16384
HOST_C0    = data2[where((data2[:,0]==1) & (data2[:,5]==1000)    ),:][0] # mpi 1 - comp 1000

MIC_C1000 = data3[where((data3[:,0]==1) & (data3[:,2]==16384) ),:][0] # mpi 1 - grid 16384
MIC_C0    = data3[where((data3[:,0]==1) & (data3[:,5]==1000)    ),:][0] # mpi 1 - comp 1000

top.set_title(str(today) + ' life_hpc2 on eurora - grid 16384')
top.grid()
top.set_xlabel('comp')
top.set_ylabel('log(time)')
top.set_yscale('log')
#top.legend()

top.plot(ACC_C1000[:,5],ACC_C1000[:,7],'-xr',HOST_C1000[:,5],HOST_C1000[:,7],'-xg',MIC_C1000[:,5],MIC_C1000[:,7],'-xc');
top.legend(('KEP20','HOST-icc','MIC'), loc = 'upper left', shadow = False, prop={'size':9})


#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc2 on eurora - COMP 1000')
bottom.grid()
bottom.set_xlabel('lattice size')
bottom.set_ylabel('time')
#bottom.set_yscale('log')
bottom.set_xlim(1000,17000)
#bottom.set_ylim(0.4,100)

bottom.plot(ACC_C0[:,3],ACC_C0[:,7],'-xr',HOST_C0[:,3],HOST_C0[:,7],'-xg',MIC_C0[:,3],MIC_C0[:,7],'-xc');
bottom.legend(('KEP20','HOST-icc','MIC'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


