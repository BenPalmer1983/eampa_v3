!#############################################################
!#  HOW TO USE:
!#
!#
!# 
!# 
!#
!#
!#
!#
!#############################################################

MODULE potential

USE kinds
USE fnc, ONLY: f, fgrad
USE cache, ONLY: clear_cache, set, get

IMPLICIT NONE

INCLUDE "potential.globals.f90"
INCLUDE "potential.interfaces.f90"

CONTAINS

INCLUDE "potential.setup.f90"
INCLUDE "potential.keys.f90"
INCLUDE "potential.search.f90"

END MODULE potential