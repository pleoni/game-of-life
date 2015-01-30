#******************************************************************************#
#* Life on EURORA                                                             *#
#******************************************************************************#

PACKAGE = life_hpc2
VERSION = 220514
RELEASE = 1

LOADPGI     = module load pgi openmpi/1.6.5--pgi--14.1
UNLOADPGI   = module unload pgi openmpi/1.6.5--pgi--14.1
LOADGNU     = module load gnu openmpi
UNLOADGNU   = module unload gnu openmpi
LOADGOMP40     = module load gomp40 openmpi-x86_64
UNLOADGOMP40   = module unload gomp40 openmpi-x86_64
LOADINTEL   = module load intel intelmpi
UNLOADINTEL = module unload intel intelmpi
UNLOADALL   = $(UNLOADPGI); $(UNLOADGNU); $(UNLOADINTEL)

TESTFILE      = test_run.dat
TESTREFERENCE = test_reference.dat
TESTPARAMS    = -r100 -c100 -s1000 -n10 -d0 -t8 # -t ignored when omp disabled

CC    = mpicc
RUN   = mpirun
KEPRUN = /opt/pgi/linux86-64/14.1/mpi/mpich/bin/mpirun
RM    = rm -f
MyO   = -O3
SIMD_PGI  = -tp=sandybridge-64 -Mvect=simd:256
SIMD_GNU  = -mavx
SIMD_INTEL = -xavx

all: pgi gnu icc mic #kep

pgi: accpgi omppgi serpgi

gnu:        ompgnu sergnu

icc:        ompicc sericc

mic:        ompmic

kep: acckep ompkep

accpgi:
	bash -c "$(UNLOADALL) ; $(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_accpgi -acc -ta=tesla:kepler -Minfo=accel -lpgacc"

omppgi: 
	bash -c "$(UNLOADALL) ; $(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omppgi $(MyO) $(SIMD_PGI) -mp=numa -fast -Minfo=vec,mp"

serpgi: 
	bash -c "$(UNLOADALL) ; $(LOADPGI) ; \
	$(CC) $(PACKAGE).c $(MyO) $(SIMD_PGI) -o $(PACKAGE)_serpgi"

ompgnu:
	bash -c "$(UNLOADALL) ; $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_ompgnu $(MyO) $(SIMD_GNU) -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1"

sergnu: 
	bash -c "$(UNLOADALL) ; $(LOADGNU) ; \
	$(CC) $(PACKAGE).c $(SIMD_GNU) $(MyO) -o $(PACKAGE)_sergnu"

ompmic:
	bash -c "$(UNLOADALL) ; $(LOADINTEL) ; \
	source $(INTEL_HOME)/bin/compilervars.sh intel64 ; \
	export I_MPI_MIC=enable ; \
	$(CC) $(MyO) $(PACKAGE).c -mmic -fopenmp -vec-report2 $(MyO) -o $(PACKAGE)_omp.mic"

ompicc:
	bash -c "$(UNLOADALL) ; $(LOADINTEL) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_ompicc $(MyO) $(SIMD_INTEL) -fopenmp -vec-report2"

omp4icc:
	bash -c "$(UNLOADALL) ; $(LOADINTEL) ; \
	source $(INTEL_HOME)/bin/compilervars.sh intel64 ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omp4icc $(MyO) $(SIMD_INTEL) -fopenmp -D OMP4 -vec-report2"

sericc: 
	bash -c "$(UNLOADALL) ; $(LOADINTEL) ; \
	$(CC) $(PACKAGE).c $(SIMD_GNU) $(MyO) $(SIMD_INTEL) -o $(PACKAGE)_sericc"

acckep:
	bash -c "$(UNLOADALL) ; $(UNLOADGOMP40) ; module load pgi ; \
	pgcc -Mmpi=mpich $(PACKAGE).c -o $(PACKAGE)_acckep -acc -ta=tesla:kepler -Minfo=accel -O3"

ompkep:
	bash -c "$(UNLOADALL) ; $(UNLOADGOMP40) ; module load pgi ; \
	pgcc -Mmpi=mpich $(PACKAGE).c -o $(PACKAGE)_ompkep $(MyO) -mp=numa -fast -Minfo=vec,mp"

omp4kep:
	bash -c "$(UNLOADALL) ; $(LOADGOMP40) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omp4kep $(MyO) -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1 -D OMP4"

clean: ; $(RM) $(PACKAGE)_accpgi  $(PACKAGE)_omppgi $(PACKAGE)_ompicc $(PACKAGE)_omp4icc $(PACKAGE)_serpgi $(PACKAGE)_ompgnu $(PACKAGE)_sergnu $(PACKAGE)_omp.mic $(PACKAGE)_acckep $(PACKAGE)_ompkep $(PACKAGE)_omp4kep $(TESTFILE) $(TESTFILE).txt

.PHONY : clean all pgi gnu mic kep

# ----------------------------------------------------------- #

tests: tests_pgi tests_gnu tests_icc

tests_pgi: test_accpgi test_omppgi test_serpgi

tests_gnu: test_ompgnu test_sergnu

tests_icc: test_ompicc test_sericc

tests_mic: test_ompmic

tests_kep: test_acckep test_ompkep

test_%gnu:  # gnu
	@bash -c "rm -f $(TESTFILE); $(UNLOADALL); $(LOADGNU); \
	$(RUN) ./$(PACKAGE)_$*gnu $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_%icc:  # intel
	@bash -c "rm -f $(TESTFILE); $(UNLOADALL); $(LOADINTEL); \
	$(RUN) ./$(PACKAGE)_$*icc $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_%kep:  # on kepler
	@bash -c "rm -f $(TESTFILE);\
	$(KEPRUN) ./$(PACKAGE)_$*kep $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_%:     # general case to match all other test targets
	@bash -c "rm -f $(TESTFILE); $(UNLOADALL); $(LOADPGI); \
	$(RUN) ./$(PACKAGE)_$* $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_reference:   # uses "sergnu" as reference, for now
	@bash -c "echo Rebuilding test reference file...; rm -i $(TESTREFERENCE); $(UNLOADALL); $(LOADGNU); \
	$(RUN) $(PACKAGE)_sergnu $(TESTPARAMS) -f$(TESTREFERENCE); echo Done."
