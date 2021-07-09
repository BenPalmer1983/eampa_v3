#!/bin/bash
python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
../kinds/src/kinds/kinds.f90 \
../sls/src/sls/sls.f90 \
../expspline/src/expspline/expspline.f90 \
src/fnc/fnc.f90 \
-m f_fnc


