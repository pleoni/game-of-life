#******************************************************************************#
#* Life on EURORA                                                             *#
#******************************************************************************#

PACKAGE = life_hpc2
VERSION = 270214
RELEASE = 1
LOADPGI = module load pgi openmpi/1.6.5--pgi--14.1
UNLOADPGI = module unload pgi openmpi/1.6.5--pgi--14.1
LOADGNU = module load gnu openmpi
UNLOADGNU = module unload gnu openmpi

CC            = mpicc
RM            = rm -fr
MyO = -O3 

all: pgiacc pgiomp gccomp

pgiacc: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_pgiacc -acc  -ta=nvidia,time -Minfo=accel -lpgacc"

pgiomp: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_pgiomp  -mp=numa -fast -mp -Minfo=vec"

pgiser: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_pgiser"

gccomp:
	bash -c "$(UNLOADPGI) ;  $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_gccomp   -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1"


clean: ; $(RM) $(PACKAGE)_pgiacc  $(PACKAGE)_pgiomp  $(PACKAGE)_gccomp  

.PHONY : clean
