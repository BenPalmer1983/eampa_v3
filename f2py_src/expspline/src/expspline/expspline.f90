!#############################################################
!#  HOW TO USE:
!#
!#
!#                                                        
!#
!#############################################################

MODULE expspline

USE kinds
USE sls, ONLY: sls_solve
USE ieee_arithmetic

IMPLICIT NONE

INCLUDE "expspline.globals.f90"
INCLUDE "expspline.interfaces.f90"

CONTAINS


INCLUDE "expspline.minverse.f90"
INCLUDE "expspline.fit.f90"


END MODULE expspline
