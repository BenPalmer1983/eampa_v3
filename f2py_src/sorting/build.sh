#!/bin/bash

python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
src/kinds/kinds.f90 \
src/bubblesort/bubblesort.f90 \
src/heapsort/heapsort.f90 \
src/mergesort/mergesort.f90 \
src/rng/rng.f90 \
src/shuffle/shuffle.f90 \
src/sort/sort.f90 \
-m f_sorting


# Test
#python3 sort_test.py
