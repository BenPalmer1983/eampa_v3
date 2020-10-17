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

MODULE efs

USE kinds
USE interp, ONLY: fill, fill_zoor, interp4dydxn
USE spline, ONLY: spline_array, spline_ab_poly, spline_ab_exp3, spline_ab_exp5

IMPLICIT NONE

INCLUDE "efs.globals.f90"
INCLUDE "efs.interfaces.f90"

CONTAINS

INCLUDE "efs.timer.f90"

INCLUDE "efs.init.f90"
INCLUDE "efs.math.f90"

INCLUDE "efs.add_config.f90"
INCLUDE "efs.nl.f90"
INCLUDE "efs.potential.f90"


INCLUDE "efs.rss_calc.f90"
INCLUDE "efs.energy.f90"
INCLUDE "efs.energy_force.f90"
INCLUDE "efs.efs.f90"

INCLUDE "efs.max_density.f90"

INCLUDE "efs.get.f90"
INCLUDE "efs.set.f90"


END MODULE efs