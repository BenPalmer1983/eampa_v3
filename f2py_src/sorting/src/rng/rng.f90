MODULE rng

USE kinds

IMPLICIT NONE




INCLUDE "rng.globals.f90"
INCLUDE "rng.interfaces.f90"

CONTAINS

INCLUDE "rng.randint.f90"
INCLUDE "rng.randint_1d.f90"
INCLUDE "rng.randint_2d.f90"
INCLUDE "rng.randdp.f90"
INCLUDE "rng.randdp_1d.f90"


END MODULE rng