MODULE bubblesort

USE kinds

IMPLICIT NONE




INCLUDE "bubblesort.globals.f90"
INCLUDE "bubblesort.interfaces.f90"

CONTAINS

INCLUDE "bubblesort.sort_1d_int.f90"
INCLUDE "bubblesort.sort_1d_dp.f90"
INCLUDE "bubblesort.sort_2d_int.f90"
INCLUDE "bubblesort.sort_2d_dp.f90"



END MODULE bubblesort