import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_acc_mpi_comm.out"
datafile2="life_host_mpi_comm.out"
datafile3="life_mic_mpi_comm.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="life_device_comm.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
#top    = fig.add_subplot(211)
bottom = fig.add_subplot(111)

#############  TOP 

ACC_COMM    = data1[where((data1[:,0]==8) & (data1[:,5]==1000) ),:][0] 
HOST_COMM    = data2[where((data2[:,0]==8) & (data2[:,5]==1000) ),:][0]
MIC_COMM    = data3[where((data3[:,0]==8) & (data3[:,5]==1000) ),:][0] 

#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc2 communication times - COMP 1000')
#bottom.grid()
bottom.set_xlabel('lattice size')
bottom.set_ylabel('time')
#bottom.set_yscale('log')
#bottom.set_xlim(1000,17000)
#bottom.set_ylim(-0.5,9)

bottom.plot(ACC_COMM[:,3],ACC_COMM[:,7],'-+r',HOST_COMM[:,3],HOST_COMM[:,7],'-+g',MIC_COMM[:,3],MIC_COMM[:,7],'-+b');
bottom.legend(('ACC-comm','HOST-comm','MIC-comm'), loc = 'upper left', shadow = False, prop={'size':8})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


