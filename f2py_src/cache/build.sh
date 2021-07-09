#!/bin/bash
python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
../kinds/src/kinds/kinds.f90 \
src/cache/cache.f90 \
-m f_cache


