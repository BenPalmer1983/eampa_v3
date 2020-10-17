#!/bin/bash
python3 -m numpy.f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" \
-lgomp \
-c \
src/kinds/kinds.f90 \
src/matrix/matrix.f90 \
src/sls/sls.f90 \
src/polyfit/polyfit.f90 \
src/interp/interp.f90 \
src/spline/spline.f90 \
src/potential/potential.f90 \
src/es/es.f90 \
-m f_es

cp *.so ../../f2py_lib/

