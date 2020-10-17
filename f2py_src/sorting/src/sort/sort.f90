MODULE sort

USE kinds

IMPLICIT NONE




INCLUDE "sort.globals.f90"
INCLUDE "sort.interfaces.f90"

CONTAINS

INCLUDE "sort.keytable.f90"
INCLUDE "sort.apply_keytable.f90"

INCLUDE "sort.sort_1d_dp_asc.f90"
INCLUDE "sort.sort_1d_dp_desc.f90"
INCLUDE "sort.sort_2d_dp_asc.f90"
INCLUDE "sort.sort_2d_dp_desc.f90"


INCLUDE "sort.sort_1d_int_asc.f90"
INCLUDE "sort.sort_1d_int_desc.f90"
INCLUDE "sort.sort_2d_int_asc.f90"
INCLUDE "sort.sort_2d_int_desc.f90"




END MODULE sort