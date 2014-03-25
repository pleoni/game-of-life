#!/bin/bash

# This script is meant to check performance scaling with the number of OpenMP threads,
# at a fixed problem size and with a fixed number of steps.
# This script can be run directly on Kepler, or sumbitted with 'qsub' on Eurora.

#PBS -N omptest
## #PBS -o omptest.out
## #PBS -e omptest.err
## #PBS -j eo
#PBS -l select=1:ncpus=16:mem=5gb:cpuspeed=2GHz
#PBS -q debug
#PBS -A CON13_INFN
## #PBS -l walltime=0:30:00

#cd /eurora/home/userexternal/ralfieri/life/openacc/game-of-life
#cd ~/game-of-life

if [[ $(hostname) != *kepler* ]]
then
  #module load pgi openmpi/1.6.5--pgi--14.1
  module load intel intelmpi mkl
  cat $PBS_NODEFILE > whereami.txt
  MYSUFFIX=ompicc
else
  MYSUFFIX=ompkep
fi

# make a temporary copy of the executable
TEMPFILE=$(mktemp)
cp ~/game-of-life/life_hpc2_$MYSUFFIX $TEMPFILE
chmod a+x $TEMPFILE

#export OMP_SCHEDULE="auto"
#echo OMP_SCHEDULE=$OMP_SCHEDULE >&2 # print to stderr
echo "OpenMP schedule: static (hardcoded)" >&2
export OMP_PROC_BIND=TRUE
echo OMP_PROC_BIND=$OMP_PROC_BIND >&2

DIM=8000
NCOMP=1000
STEPS=10
echo DIM=$DIM, NCOMP=$NCOMP, STEPS=$STEPS >&2

for T in {1..16}
do
  CMD="$TEMPFILE -r$DIM -c$DIM -n$NCOMP -s$STEPS -d0 -t$T -a"
  echo "# $CMD"
  eval $CMD
done

DIM=8000
NCOMP=0
STEPS=100
echo DIM=$DIM, NCOMP=$NCOMP, STEPS=$STEPS >&2

for T in {1..16}
do
  CMD="$TEMPFILE -r$DIM -c$DIM -n$NCOMP -s$STEPS -d0 -t$T -a"
  echo "# $CMD"
  eval $CMD
done

DIM=4000
NCOMP=0
STEPS=1000
echo DIM=$DIM, NCOMP=$NCOMP, STEPS=$STEPS >&2

for T in {1..16}
do
  CMD="$TEMPFILE -r$DIM -c$DIM -n$NCOMP -s$STEPS -d0 -t$T -a"
  echo "# $CMD"
  eval $CMD
done

# delete the temporary file
rm $TEMPFILE
