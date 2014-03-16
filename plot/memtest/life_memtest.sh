#!/bin/bash

# These timings should give an idea of how much influence the way we allocate memory
# has on memory copy times and overall program execution time.
# This script can be run directly on Kepler, or sumbitted with 'qsub' on Eurora.

#PBS -N life_memtest
#PBS -o life_memtest.out
#PBS -e life_memtest.err
#PBS -l select=1:ncpus=8:ngpus=2:mem=5gb
#PBS -q debug
#PBS -A CON13_INFN
##PBS -l walltime=0:15:00

#cd /eurora/home/userexternal/ralfieri/life/openacc/game-of-life
cd ~/game-of-life

if [[ $(hostname) != *kepler* ]]
then
  module load pgi openmpi/1.6.5--pgi--14.1
  echo $PBS_NODEFILE > whereami.txt
  MYSUFFIX=accpgi
else
  MYSUFFIX=acckep
fi


# Allocazione riga per riga
for DIM in $(echo 2000 4000 6000 8000 10000 12000 14000 16000)
do
  CMD="./life_hpc2_$MYSUFFIX -r$DIM -c$DIM -n0 -s5 -d0"
  echo "# $CMD"
  eval $CMD
done

echo ""

# Allocazione contigua
for DIM in $(echo 2000 4000 6000 8000 10000 12000 14000 16000)
do
  CMD="./life_hpc2_$MYSUFFIX -r$DIM -c$DIM -n0 -s5 -d0 -a"
  echo "# $CMD"
  eval $CMD
done
