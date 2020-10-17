#!/bin/bash
python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
src/kinds/kinds.f90 \
src/interp/interp.f90 \
-m f_interp

# Test
python3 interp.py
