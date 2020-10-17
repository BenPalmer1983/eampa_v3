#!/bin/bash
#f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -lmatrix -c kinds.f90 sls.f90 -m f_sls

f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -c kinds.f90 sls.f90 -m f_sls f_matrix.so




f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -lf_matrix -c kinds.f90 sls.f90 -m f_sls

f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -lf_matrix -c kinds.f90 sls.f90 -m f_sls

f2py -m f_sls -h sgnFile.pyf sls.f90



f2py -m multMatrices -h sgnFile.pyf matrixMult.F90
f2py -c --fcompiler=gnu95 sgnFile.pyf matrixMult.F90 -L{Full_Path_to_Library} -lmatMultiplication

f2py -m f_sls -h sgnFile.pyf sls.f90


-L{Full_Path_to_Library} -lmatMultiplication


 
f2py -c --fcompiler=gnu95 kinds.f90 sls.f90  -L"/cloud/Code/lib" -lf_matrix -m f_sls 



f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -lf_matrix -c kinds.f90 sls.f90 -m f_sls



f2py --quiet --f90flags="-fopenmp -O3 -ffast-math" -L/cloud/Code/lib -lgomp -lmatrix -c kinds.f90 sls.f90 -m f_sls


gfortran -O3 kinds.f90 sls.f90 -L/cloud/Code/lib -lgomp -lmatrix -o sls.o 