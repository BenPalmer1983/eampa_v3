MODULE matrix

USE kinds

IMPLICIT NONE

INCLUDE "matrix.interfaces.f90"
INCLUDE "matrix.globals.f90"

CONTAINS

INCLUDE "matrix.pivot.f90"
INCLUDE "matrix.inv.f90"



END MODULE matrix
