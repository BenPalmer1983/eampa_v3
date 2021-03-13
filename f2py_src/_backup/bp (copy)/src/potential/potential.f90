!#############################################################
!#  HOW TO USE:
!#
!#
!# interp.interp4(x, x_points, y_points, y)
!# interp.interp4dydx(x, x_points, y_points, y)
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

IMPLICIT NONE

INCLUDE "potential.globals.f90"
INCLUDE "potential.interfaces.f90"

CONTAINS

INCLUDE "potential.add.f90"

END MODULE potential