#!/bin/bash

# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=5
export PROC_COUNT=1


module purge; module load bluebear
module load bear-apps/2018a
module load imkl 2019.5.281-gompi-2019b
module load Python 3.7.4-GCCcore-8.3.0
module load matplotlib 3.1.1-foss-2019b-Python-3.7.4


thisdir=$(pwd)
cd ../../
topdir=$(pwd)
./package.sh
cd $thisdir
python3 $topdir/eampa.py input.in

