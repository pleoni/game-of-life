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

all: accpgi omppgi ompgnu

accpgi:
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_accpgi -acc  -ta=nvidia,time -Minfo=accel -lpgacc"

omppgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omppgi  -mp=numa -fast -mp -Minfo=vec"

serpgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_serpgi"

ompgnu:
	bash -c "$(UNLOADPGI) ;  $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_ompgnu   -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1"


clean: ; $(RM) $(PACKAGE)_accpgi  $(PACKAGE)_omppgi  $(PACKAGE)_ompgnu 

.PHONY : clean
