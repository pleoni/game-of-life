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

plotfile="plot_mpi_comp_1000.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(311)
middle = fig.add_subplot(312)
bottom = fig.add_subplot(313)

M8ACC = data1[where((data1[:,0]==8) & (data1[:,5]==1000) ),:][0] # mpi 8 - Comp 1000
M4ACC = data1[where((data1[:,0]==4) & (data1[:,5]==1000) ),:][0] # mpi 4 - Comp 1000
M2ACC = data1[where((data1[:,0]==2) & (data1[:,5]==1000) ),:][0] # mpi 2 - comp 1000
M1ACC = data1[where((data1[:,0]==1) & (data1[:,5]==1000) ),:][0] # mpi 1 - Comp 1000

M8HOST = data2[where((data2[:,0]==8) & (data2[:,5]==1000) ),:][0] # mpi 8 - Comp 1000
M4HOST = data2[where((data2[:,0]==4) & (data2[:,5]==1000) ),:][0] # mpi 4 - Comp 1000
M2HOST = data2[where((data2[:,0]==2) & (data2[:,5]==1000) ),:][0] # mpi 2 - comp 1000
M1HOST = data2[where((data2[:,0]==1) & (data2[:,5]==1000) ),:][0] # mpi 1 - Comp 1000

M8MIC = data3[where((data3[:,0]==8) & (data3[:,5]==1000) ),:][0] # mpi 8 - Comp 1000
M4MIC = data3[where((data3[:,0]==4) & (data3[:,5]==1000) ),:][0] # mpi 4 - Comp 1000
M2MIC = data3[where((data3[:,0]==2) & (data3[:,5]==1000) ),:][0] # mpi 2 - comp 1000
M1MIC = data3[where((data3[:,0]==1) & (data3[:,5]==1000) ),:][0] # mpi 1 - Comp 1000


#############  TOP 

top.set_title(str(today) + ' LIFE MPI - K20 on eurora - NCOMP=1000')
top.grid()
top.set_xlabel('lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

top.plot(M8ACC[:,3],M8ACC[:,8],'-xr',M4ACC[:,3],M4ACC[:,8],'-xg',M1ACC[:,3],M1ACC[:,8],'-xc');
top.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})

#############  MIDDLE

middle.set_title(str(today) + ' LIFE MPI - HOST on eurora - NCOMP=1000')
middle.grid()
middle.set_xlabel('lattice Size')
middle.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

middle.plot(M8HOST[:,3],M8HOST[:,8],'-xr',M4HOST[:,3],M4HOST[:,8],'-xg',M1HOST[:,3],M1HOST[:,8],'-xc');
middle.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})

#############  BOTTOM 


bottom.set_title(str(today) + ' LIFE MPI - MIC on eurora - NCOMP=1000')
bottom.grid()
bottom.set_xlabel('lattice size')
bottom.set_ylabel('time')

bottom.plot(M8MIC[:,3],M8MIC[:,8],'-xr',M4MIC[:,3],M4MIC[:,8],'-xg',M1MIC[:,3],M1MIC[:,8],'-xc');
bottom.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.8)

plt.savefig(plotfile)
#plt.show()


