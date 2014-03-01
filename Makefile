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
TESTPARAMS    = -r100 -c100 -s1000 -n10

CC    = mpicc
RUN   = mpirun
RM    = rm -f
MyO   = -O3

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
	$(CC) $(MyO) $(PACKAGE).c -mmic -fopenmp -vec-report1 -o $(PACKAGE)_omp.mic"

clean: ; $(RM) $(PACKAGE)_accpgi  $(PACKAGE)_omppgi  $(PACKAGE)_ompgnu $(PACKAGE)_omp.mic $(TESTFILE)

.PHONY : clean


tests: test_accpgi test_omppgi test_ompgnu

test_accpgi:  $(PACKAGE)_accpgi
	@bash -c "rm -f $(TESTFILE); $(LOADPGI); \
	$(RUN) $(PACKAGE)_accpgi $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi"

test_omppgi:  $(PACKAGE)_omppgi
	@bash -c "rm -f $(TESTFILE);  $(LOADPGI); \
	$(RUN) $(PACKAGE)_ompgnu $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi"

test_ompgnu:  $(PACKAGE)_ompgnu
	@bash -c "rm -f $(TESTFILE);  $(UNLOADPGI); $(LOADGNU); \
	$(RUN) $(PACKAGE)_ompgnu $(TESTPARAMS) -f$(TESTFILE); \
	echo Comparing $(TESTFILE) to $(TESTREFERENCE)...; \
	if diff $(TESTFILE) $(TESTREFERENCE) &> /dev/null; then echo --- $@ OK ---; else echo --- $@ FAIL ---; fi"

test_reference:   # uses ompgnu as reference, for now
	@bash -c "echo Rebuilding test reference file...; rm -i $(TESTREFERENCE); $(LOADGNU); \
	$(RUN) $(PACKAGE)_ompgnu $(TESTPARAMS) -f$(TESTREFERENCE); echo Done."
