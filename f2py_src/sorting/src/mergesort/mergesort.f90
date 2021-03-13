!############################################
!#  
!#  Example usage:
!#  
!#  arr = rng.randint_1d(1,50,1000000)
!#  mergesort.set_1d_int(arr)
!#  mergesort.sort_1d_int()
!#  arr2 = mergesort.arr_int_1d
!#  kt = mergesort.key_table
!#  
!############################################


MODULE mergesort

USE kinds

IMPLICIT NONE




INCLUDE "mergesort.globals.f90"
INCLUDE "mergesort.interfaces.f90"

CONTAINS

INCLUDE "mergesort.set.f90"
INCLUDE "mergesort.return.f90"
INCLUDE "mergesort.sort.f90"
INCLUDE "mergesort.sort_1d_int.f90"
INCLUDE "mergesort.sort_2d_int.f90"
INCLUDE "mergesort.sort_1d_dp.f90"
INCLUDE "mergesort.sort_2d_dp.f90"

INCLUDE "mergesort.match_key_table.f90"


END MODULE mergesort