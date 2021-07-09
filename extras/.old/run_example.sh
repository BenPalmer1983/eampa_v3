#!/bin/bash
export PYTHONPATH=$PYTHONPATH:"/home/ben/pylib"

# Set the number of threads to 1
export OMP_NUM_THREADS=5

# Run
python3 example_spline.py

