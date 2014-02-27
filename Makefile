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
LOADINTEL = module load intel intelmpi
UNLOADINTEL = module unload intel intelmpi

CC            = mpicc
RM            = rm -fr
MyO = -O3 

all: accpgi omppgi ompgnu

accpgi:
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_accpgi -acc  -ta=nvidia -Minfo=accel -lpgacc"

omppgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omppgi  -mp=numa -fast -mp -Minfo=vec"

serpgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_serpgi"

ompgnu:
	bash -c "$(UNLOADPGI) ;  $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_ompgnu   -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1"

ompmic:
	bash -c "$(UNLOADGNU) ;  $(LOADINTEL) ; \
	source $(INTEL_HOME)/bin/compilervars.sh intel64 ; \
	export I_MPI_MIC=enable ; \
	export LD_LIBRARY_PATH= ; \
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/eurora/prod/compilers/intel/cs-xe-2013/binary/composer_xe_2013/lib/mic ; \
	mpicc $(MyO) $(PACKAGE).c -mmic -fopenmp -vec-report1 -o $(PACKAGE)_omp.mic"

clean: ; $(RM) $(PACKAGE)_accpgi  $(PACKAGE)_omppgi  $(PACKAGE)_ompgnu $(PACKAGE)_omp.mic

.PHONY : clean
