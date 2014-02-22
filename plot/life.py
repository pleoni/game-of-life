import matplotlib.pyplot as plt
from numpy import *

datafile="life.out"
plotfile=datafile+".png"
data = loadtxt(datafile)

fig = plt.figure() # apre una nuova figura
top    = fig.add_subplot(211)
bottom = fig.add_subplot(212)

#############  TOP 

M8C1000 = data[where((data[:,0]==8) & (data[:,5]==1000) ),:][0] #
M8C100  = data[where((data[:,0]==8) & (data[:,5]==100)  ),:][0] # 
M8C10   = data[where((data[:,0]==8) & (data[:,5]==10)   ),:][0] # 
M8C0    = data[where((data[:,0]==8) & (data[:,5]==0)    ),:][0] # 

M4C1000 = data[where((data[:,0]==4) & (data[:,5]==1000) ),:][0] #
M4C100  = data[where((data[:,0]==4) & (data[:,5]==100)  ),:][0] # 
M4C10   = data[where((data[:,0]==4) & (data[:,5]==10)   ),:][0] # 
M4C0    = data[where((data[:,0]==4) & (data[:,5]==0)    ),:][0] # 

M2C1000 = data[where((data[:,0]==2) & (data[:,5]==1000) ),:][0] #
M2C100  = data[where((data[:,0]==2) & (data[:,5]==100)  ),:][0] # 
M2C10   = data[where((data[:,0]==2) & (data[:,5]==10)   ),:][0] # 
M2C0    = data[where((data[:,0]==2) & (data[:,5]==0)    ),:][0] # 

M1C1000 = data[where((data[:,0]==1) & (data[:,5]==1000) ),:][0] #
M1C100  = data[where((data[:,0]==1) & (data[:,5]==100)  ),:][0] # 
M1C10   = data[where((data[:,0]==1) & (data[:,5]==10)   ),:][0] # 
M1C0    = data[where((data[:,0]==1) & (data[:,5]==0)    ),:][0] # 

top.set_title('life_hpc_acc_mpi (eurora, feb2014) - NCOMP=1000 ')
top.grid()
top.set_xlabel('Lattice Size')
top.set_ylabel('time')
#top.set_yscale('log') 
#top.legend()

top.plot(M8C1000[:,3],M8C1000[:,8],'-xr',M4C1000[:,3],M4C1000[:,8],'-xg',M1C1000[:,3],M1C1000[:,8],'-xc');
top.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})


#############  BOTTOM 


bottom.set_title('NCOMP=0')
bottom.grid()
bottom.set_xlabel('Lattice size')
bottom.set_ylabel('time')

bottom.plot(M8C0[:,3],M8C0[:,8],'-xr',M4C0[:,3],M4C0[:,8],'-xg',M1C0[:,3],M1C0[:,8],'-xc');
bottom.legend(('MPI-8','MPI-4','MPI-1'), loc = 'upper left', shadow = False, prop={'size':9})


plt.subplots_adjust(hspace=0.5)

plt.savefig(plotfile)
#plt.show()


