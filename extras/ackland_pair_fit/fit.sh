#!/bin/bash
export PYTHONPATH=$PYTHONPATH:"/home/ben/pylib"

# Set the number of threads to 1
export OMP_NUM_THREADS=4

# Run
python3 ackland_pair_fit.py "$@"

