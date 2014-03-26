import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys
import datetime

datafile="life_mic_mpi.out"
if len(sys.argv) > 1:
    datafile=sys.argv[1]

plotfile=datafile+".png"
data = loadtxt(datafile)

today = datetime.date.today()

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(211)
bottom = fig.add_subplot(212)

#############  TOP 

M8C1000 = data[where((data[:,0]==8) & (data[:,5]==1000) ),:][0] # mpi 8 - Comp 1000
M8C100  = data[where((data[:,0]==8) & (data[:,5]==100)  ),:][0] # mpi 8 - Comp 100
M8C10   = data[where((data[:,0]==8) & (data[:,5]==10)   ),:][0] # mpi 8 - Comp 10
M8C0    = data[where((data[:,0]==8) & (data[:,5]==0)    ),:][0] # mpi 8 - Comp 0

M4C1000 = data[where((data[:,0]==4) & (data[:,5]==1000) ),:][0] # mpi 4 - Comp 1000
M4C100  = data[where((data[:,0]==4) & (data[:,5]==100)  ),:][0] # mpi 4 - comp 100
M4C10   = data[where((data[:,0]==4) & (data[:,5]==10)   ),:][0] # mpi 4 - Comp 10
M4C0    = data[where((data[:,0]==4) & (data[:,5]==0)    ),:][0] # mpi 4 - comp 0 

M2C1000 = data[where((data[:,0]==2) & (data[:,5]==1000) ),:][0] # mpi 2 - comp 1000
M2C100  = data[where((data[:,0]==2) & (data[:,5]==100)  ),:][0] # mpi 2 - comp 100
M2C10   = data[where((data[:,0]==2) & (data[:,5]==10)   ),:][0] # mpi 2 - Comp 10
M2C0    = data[where((data[:,0]==2) & (data[:,5]==0)    ),:][0] # mpi 2 - comp 0

M1C1000 = data[where((data[:,0]==1) & (data[:,5]==1000) ),:][0] # mpi 1 - Comp 1000
M1C100  = data[where((data[:,0]==1) & (data[:,5]==100)  ),:][0] # mpi 1 - Comp 100
M1C10   = data[where((data[:,0]==1) & (data[:,5]==10)   ),:][0] # mpi 1 - Comp 10
M1C0    = data[where((data[:,0]==1) & (data[:,5]==0)    ),:][0] # mpi 1 - comp 0

top.set_title(str(today) + ' life_hpc_mic_mpi on eurora - NCOMP=1000')
top.grid()
top.set_xlabel('Lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

top.plot(M8C1000[:,3],M8C1000[:,8],'-xr',M4C1000[:,3],M4C1000[:,8],'-xg',M1C1000[:,3],M1C1000[:,8],'-xc');
top.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})


#############  BOTTOM 


bottom.set_title(str(today) + ' life_hpc_mic_mpi on eurora - NCOMP=0')
bottom.grid()
bottom.set_xlabel('Lattice size')
bottom.set_ylabel('time')

bottom.plot(M8C0[:,3],M8C0[:,8],'-xr',M4C0[:,3],M4C0[:,8],'-xg',M1C0[:,3],M1C0[:,8],'-xc');
bottom.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


