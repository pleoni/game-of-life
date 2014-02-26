#******************************************************************************#
#* Life_hpc	                                                              *#
#******************************************************************************#
PACKAGE = life_hpc2
VERSION = 260214
RELEASE = 1

CC            = mpicc
RM            = rm -fr
MyO = -O3 

all: acc omp serial


acc: 
	$(CC) $(PACKAGE).c -o $(PACKAGE)_acc.exe -acc -DCOMP -ta=nvidia,time -Minfo=accel -lpgacc

omp: 
	$(CC) $(PACKAGE).c -o $(PACKAGE)_omp.exe -DCOMP -DOMP -mp=numa -fast -mp -Minfo=vec 

serial: 
	$(CC) $(PACKAGE).c -o $(PACKAGE)_serial.exe -DCOMP

clean: ; $(RM) *.exe

.PHONY : clean
