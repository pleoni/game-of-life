import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_acc_mpi.out"
datafile2="life_host_mpi.out"
datafile3="life_mic_mpi.out"
datafile4="life_acc_mpi_comm.out"
datafile5="life_host_mpi_comm.out"
datafile6="life_mic_mpi_comm.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="plot_eurora_comp_comm.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)
data4 = loadtxt(datafile4)
data5 = loadtxt(datafile5)
data6 = loadtxt(datafile6)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(311)
middle = fig.add_subplot(312)
bottom = fig.add_subplot(313)

ACC_COMP    = data1[where((data1[:,0]==8) & (data1[:,5]==1000) ),:][0] 
HOST_COMP    = data2[where((data2[:,0]==8) & (data2[:,5]==1000) ),:][0] 
MIC_COMP    = data3[where((data3[:,0]==8) & (data3[:,5]==1000) ),:][0] 
ACC_COMM    = data4[where((data4[:,0]==8) & (data4[:,5]==1000) ),:][0] 
HOST_COMM    = data5[where((data5[:,0]==8) & (data5[:,5]==1000) ),:][0]
MIC_COMM    = data6[where((data6[:,0]==8) & (data6[:,5]==1000) ),:][0] 

#############  TOP 

top.set_title(str(today) + ' LIFE - Comp vs Comm - MPI 8 - K20')
top.grid()
top.set_xlabel('lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()
#top.set_xlim(-10,1050)
top.plot(ACC_COMP[:,3],ACC_COMP[:,7],'-xr',ACC_COMM[:,3],ACC_COMM[:,7],':+r');
top.grid(False)

top.legend(('ACC','ACC-comm'), loc = 'upper left', shadow = False, prop={'size':8})


#############  MIDDLE

middle.set_title(str(today) + ' LIFE - Comp vs Comm - MPI 8 - HOST')
middle.grid()
middle.set_xlabel('lattice Size')
middle.set_ylabel('time')
#middle.set_yscale('log') 
#middle.set_xlim(-10,1050)
middle.set_ylim(-20,250)
middle.grid(False)

middle.plot(HOST_COMP[:,3],HOST_COMP[:,7],'-xg',HOST_COMM[:,3],HOST_COMM[:,7],':+g');
middle.legend(('HOST','HOST-comm'), loc = 'upper left', shadow = False, prop={'size':8})


#############  BOTTOM

bottom.set_title(str(today) + ' LIFE - Comp vs Comm - MPI 8 - MIC')
bottom.grid()
bottom.set_xlabel('lattice Size')
bottom.set_ylabel('time')
#bottom.set_yscale('log') 
bottom.set_ylim(-20,180)
bottom.grid(False)

bottom.plot(MIC_COMP[:,3],MIC_COMP[:,7],'-xb',MIC_COMM[:,3],MIC_COMM[:,7],':+b');
bottom.legend(('MIC','MIC-comm'), loc = 'upper left', shadow = False, prop={'size':8})

plt.subplots_adjust(hspace=0.8)

plt.savefig(plotfile)
#plt.show()


