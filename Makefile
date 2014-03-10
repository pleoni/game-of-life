#******************************************************************************#
#* Life on EURORA                                                             *#
#******************************************************************************#

PACKAGE = life_hpc2
VERSION = 270214
RELEASE = 1

LOADPGI     = module load pgi openmpi/1.6.5--pgi--14.1
UNLOADPGI   = module unload pgi openmpi/1.6.5--pgi--14.1
LOADGNU     = module load gnu openmpi
UNLOADGNU   = module unload gnu openmpi
LOADINTEL   = module load intel intelmpi
UNLOADINTEL = module unload intel intelmpi
UNLOADALL   = $(UNLOADPGI); $(UNLOADGNU); $(UNLOADINTEL)

TESTFILE      = test_run.dat
TESTREFERENCE = test_reference.dat
TESTPARAMS    = -r100 -c100 -s1000 -n10 -d0

CC    = mpicc
RUN   = mpirun
KEPRUN = /opt/pgi/linux86-64/14.1/mpi/mpich/bin/mpirun
RM    = rm -f
MyO   = -O3

all: pgi gnu mic kep

pgi: accpgi omppgi serpgi

gnu:        ompgnu sergnu

mic:        ompmic

kep: acckep ompkep

accpgi:
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_accpgi -acc -ta=nvidia -Minfo=accel -lpgacc"

omppgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omppgi $(MyO) -mp=numa -fast -mp -Minfo=vec"

serpgi: 
	bash -c "$(LOADPGI) ; \
	$(CC) $(PACKAGE).c $(MyO) -o $(PACKAGE)_serpgi"

ompgnu:
	bash -c "$(UNLOADPGI) ;  $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_ompgnu $(MyO) -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=1"

sergnu: 
	bash -c "$(UNLOADPGI) ;  $(LOADGNU) ; \
	$(CC) $(PACKAGE).c -o $(PACKAGE)_sergnu"

ompmic:
	bash -c "$(UNLOADGNU) ;  $(LOADINTEL) ; \
	source $(INTEL_HOME)/bin/compilervars.sh intel64 ; \
	export I_MPI_MIC=enable ; \
	$(CC) $(MyO) $(PACKAGE).c -mmic -fopenmp -vec-report2 -o $(PACKAGE)_omp.mic"

acckep:
	pgcc -Mmpi=mpich $(PACKAGE).c -o $(PACKAGE)_acckep -acc -ta=nvidia -Minfo=accel -lpgacc

ompkep:
	pgcc -Mmpi=mpich $(PACKAGE).c -o $(PACKAGE)_ompkep $(MyO) -mp=numa -fast -mp -Minfo=vec


clean: ; $(RM) $(PACKAGE)_accpgi  $(PACKAGE)_omppgi $(PACKAGE)_serpgi $(PACKAGE)_ompgnu $(PACKAGE)_sergnu $(PACKAGE)_omp.mic $(PACKAGE)_acckep $(PACKAGE)_ompkep $(TESTFILE)

.PHONY : clean all pgi gnu mic kep

# ----------------------------------------------------------- #

tests: tests_pgi tests_gnu

tests_pgi: test_accpgi test_omppgi test_serpgi

tests_gnu: test_ompgnu test_sergnu

tests_mic: test_ompmic

tests_kep: test_acckep test_ompkep

test_%kep:  # special case to match on kepler
	@bash -c "rm -f $(TESTFILE);\
	$(KEPRUN) ./$(PACKAGE)_$*kep $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_%:     # general case to match all other test targets
	@bash -c "rm -f $(TESTFILE); $(LOADPGI); \
	$(RUN) ./$(PACKAGE)_$* $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi; echo "

test_reference:   # uses "sergnu" as reference, for now
	@bash -c "echo Rebuilding test reference file...; rm -i $(TESTREFERENCE); $(UNLOADALL); $(LOADGNU); \
	$(RUN) $(PACKAGE)_sergnu $(TESTPARAMS) -f$(TESTREFERENCE); echo Done."
