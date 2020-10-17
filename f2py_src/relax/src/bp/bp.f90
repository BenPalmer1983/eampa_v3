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

MODULE bp

USE kinds
USE matrix, ONLY: inv
USE polyfit, ONLY: fit
USE interp, ONLY: fill, fill_zoor, interp4dydxn
USE spline, ONLY: spline_array, spline_ab_poly, spline_ab_exp3, spline_ab_exp5
USE sls, ONLY: solve

IMPLICIT NONE

INCLUDE "bp.globals.f90"
INCLUDE "bp.interfaces.f90"

CONTAINS

INCLUDE "bp.init.f90"
INCLUDE "bp.math.f90"
INCLUDE "bp.bm_fit.f90"


INCLUDE "bp.add_config.f90"
INCLUDE "bp.add_calculation.f90"
INCLUDE "bp.add_known.f90"
INCLUDE "bp.nl.f90"
INCLUDE "bp.potential.f90"

INCLUDE "bp.energy.f90"

INCLUDE "bp.bp.f90"

INCLUDE "bp.set_rss.f90"

!INCLUDE "efs.get.f90"



END MODULE bp