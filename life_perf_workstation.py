### Dependencies: yum -y install python-matplotlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile="my_host.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="device_perf_workstation.png"
data = loadtxt(datafile)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(311)
middle = fig.add_subplot(312)
bottom = fig.add_subplot(313)

#############  TOP 

HOST_GRID = data[where((data[:,1]==4) & (data[:,2]==8192) ),:][0] # mpi 1 - grid 16384
HOST_COMP = data[where((data[:,1]==4) & (data[:,5]==500)    ),:][0] # mpi 1 - comp 1000
HOST_THREAD = data[where((data[:,2]==8192) & (data[:,5]==500)    ),:][0] # mpi 1 - comp 1000

top.set_title(str(today) + ' life_hpc2 on workstation - GRID 8192')
top.grid()
top.set_xlabel('Comp')
top.set_ylabel('Time')
#top.set_yscale('log') 
#top.legend()

top.plot(HOST_GRID[:,5],HOST_GRID[:,7],'-xr');
#top.legend(('ciao'), loc = 'upper left', shadow = False, prop={'size':9})

#############  MIDDLE


middle.set_title(str(today) + ' life_hpc2 on workstation - COMP 500')
middle.grid()
middle.set_xlabel('Lattice Size')
middle.set_ylabel('Time')
#middle.set_yscale('log')
#middle.set_xlim(1000,17000)
#middle.set_ylim(0.4,100)

middle.plot(HOST_COMP[:,3],HOST_COMP[:,7],'-xr');
#middle.legend(('localhost'), loc = 'upper left', shadow = False, prop={'size':9})

#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc2 on workstation - GRID 8192 - THREADS')
bottom.grid()
bottom.set_xlabel('Threads')
bottom.set_ylabel('Time')
#bottom.set_yscale('log')
#bottom.set_xlim(1000,17000)
#bottom.set_ylim(0.4,100)
bottom.plot(HOST_THREAD[:,1],HOST_THREAD[:,7],'-xr');
#bottom.legend(('localhost'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.8)

plt.savefig(plotfile)
#plt.show()


