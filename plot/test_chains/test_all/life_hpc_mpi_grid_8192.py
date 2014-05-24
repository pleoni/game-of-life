import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_acc_mpi.out"
datafile2="life_host_mpi.out"
datafile3="life_mic_mpi.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="plot_mpi_test.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(311)
middle = fig.add_subplot(312)
bottom = fig.add_subplot(313)

#############  TOP 

M8ACC    = data1[where((data1[:,0]==8) & (data1[:,2]==8192)    ),:][0] # mpi 8 - grid 8192
M4ACC    = data1[where((data1[:,0]==4) & (data1[:,2]==8192)    ),:][0] # mpi 4 - grid 8192
M2ACC    = data1[where((data1[:,0]==2) & (data1[:,2]==8192)    ),:][0] # mpi 2 - grid 8192
M1ACC    = data1[where((data1[:,0]==1) & (data1[:,2]==8192)    ),:][0] # mpi 1 - grid 8192

top.set_title(str(today) + ' LIFE - MPI on eurora - K20 - grid 8192')
top.grid()
top.set_xlabel('Comp')
top.set_ylabel('Time(log)')
#top.set_yscale('log') 
#top.legend()
top.set_xlim(-10,1050)
top.plot(M8ACC[:,5],M8ACC[:,8],'-xr',M4ACC[:,5],M4ACC[:,8],'-+g',M1ACC[:,5],M1ACC[:,8],'-oc');

top.legend(('MPI-8','MPI-4','MPI-1'), loc = 'lower right', shadow = False, prop={'size':8})


#############  MIDDLE
M8HOST    = data2[where((data2[:,0]==8) & (data2[:,2]==8192)    ),:][0] # mpi 8 - grid 8192
M4HOST    = data2[where((data2[:,0]==4) & (data2[:,2]==8192)    ),:][0] # mpi 4 - grid 8192
M2HOST    = data2[where((data2[:,0]==2) & (data2[:,2]==8192)    ),:][0] # mpi 2 - grid 8192
M1HOST    = data2[where((data2[:,0]==1) & (data2[:,2]==8192)    ),:][0] # mpi 1 - grid 8192

middle.set_title(str(today) + ' LIFE - MPI on eurora - HOST - grid 8192')
middle.grid()
middle.set_xlabel('Comp')
middle.set_ylabel('Time(log)')
#middle.set_yscale('log') 
middle.set_xlim(-10,1050)

middle.plot(M8HOST[:,5],M8HOST[:,8],'-xr',M4HOST[:,5],M4HOST[:,8],'-+g',M1HOST[:,5],M1HOST[:,8],'-oc');
middle.legend(('MPI-8','MPI-4','MPI-1'), loc = 'lower right', shadow = False, prop={'size':8})


#############  BOTTOM
M8MIC    = data3[where((data3[:,0]==8) & (data3[:,2]==8192)    ),:][0] # mpi 8 - grid 8192
M4MIC    = data3[where((data3[:,0]==4) & (data3[:,2]==8192)    ),:][0] # mpi 4 - grid 8192
M2MIC    = data3[where((data3[:,0]==2) & (data3[:,2]==8192)    ),:][0] # mpi 2 - grid 8192
M1MIC    = data3[where((data3[:,0]==1) & (data3[:,2]==8192)    ),:][0] # mpi 1 - grid 8192
bottom.set_title(str(today) + 'LIFE - MPI on eurora - MIC - grid 8192')
bottom.grid()
bottom.set_xlabel('Comp')
bottom.set_ylabel('Time(log)')
#bottom.set_yscale('log') 
bottom.set_xlim(-10,1050)

bottom.plot(M8MIC[:,5],M8MIC[:,8],'-xr',M4MIC[:,5],M4MIC[:,8],'-+g',M1MIC[:,5],M1MIC[:,8],'-oc');
bottom.legend(('MPI-8','MPI-4','MPI-1'), loc = 'lower right', shadow = False, prop={'size':8})

plt.subplots_adjust(hspace=0.8)

plt.savefig(plotfile)
#plt.show()


