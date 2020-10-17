MODULE shuffle

USE kinds
USE rng, ONLY: randint

IMPLICIT NONE




INCLUDE "shuffle.globals.f90"
INCLUDE "shuffle.interfaces.f90"

CONTAINS

INCLUDE "shuffle.shuffle_1d_int.f90"


END MODULE shuffle