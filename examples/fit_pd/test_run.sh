#!/bin/bash

# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=4
export PROC_COUNT=1

python3 ../../eampa.py input.in

