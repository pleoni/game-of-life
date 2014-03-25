import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile1="life_host_icc.out"
datafile2="life_host_gnu.out"
datafile3="life_host_pgi.out"

if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile="compilers_perf_eurora.png"
data1 = loadtxt(datafile1)
data2 = loadtxt(datafile2)
data3 = loadtxt(datafile3)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(211)
bottom = fig.add_subplot(212)

#############  TOP 

ICC_C1000 = data1[where((data1[:,0]==1) & (data1[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
ICC_C0    = data1[where((data1[:,0]==1) & (data1[:,5]==0)    ),:][0] # mpi 1 - comp 0

GNU_C1000 = data2[where((data2[:,0]==1) & (data2[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
GNU_C0    = data2[where((data2[:,0]==1) & (data2[:,5]==0)    ),:][0] # mpi 1 - comp 0

PGI_C1000 = data3[where((data3[:,0]==1) & (data3[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
PGI_C0    = data3[where((data3[:,0]==1) & (data3[:,5]==0)    ),:][0] # mpi 1 - comp 0

top.set_title(str(today) + ' life_hpc2 on eurora - NCOMP=1000')
top.grid()
top.set_xlabel('Lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

top.plot(ICC_C1000[:,3],ICC_C1000[:,8],'-xr',GNU_C1000[:,3],GNU_C1000[:,8],'-xg',PGI_C1000[:,3],PGI_C1000[:,8],'-xc');
top.legend(('icc','gnu','pgi'), loc = 'upper left', shadow = False, prop={'size':9})


#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc2 on eurora - NCOMP=0')
bottom.grid()
bottom.set_xlabel('Lattice size')
bottom.set_ylabel('time')

bottom.plot(ICC_C0[:,3],ICC_C0[:,8],'-xr',GNU_C0[:,3],GNU_C0[:,8],'-xg',PGI_C0[:,3],PGI_C0[:,8],'-xc');
bottom.legend(('icc','gnu','pgi'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


