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
!#
!#
!#############################################################

MODULE es

USE kinds
USE matrix, ONLY: inv
USE polyfit, ONLY: fit
USE interp, ONLY: fill, fill_zoor, interp4dydxn
USE spline, ONLY: spline_array, spline_ab_poly, spline_ab_exp3, spline_ab_exp5
USE sls, ONLY: solve

IMPLICIT NONE

INCLUDE "es.globals.f90"
INCLUDE "es.interfaces.f90"

CONTAINS

INCLUDE "es.init.f90"
INCLUDE "es.math.f90"

INCLUDE "es.potential.f90"
INCLUDE "es.add_config.f90"
INCLUDE "es.add_calculation.f90"


INCLUDE "es.nl.f90"
INCLUDE "es.energy.f90"
INCLUDE "es.force.f90"


INCLUDE "es.es.f90"


!INCLUDE "bp.add_config.f90"
!INCLUDE "bp.add_calculation.f90"
!INCLUDE "bp.add_known.f90"
!INCLUDE "bp.nl.f90"

!INCLUDE "bp.energy.f90"

!INCLUDE "bp.bp.f90"

!INCLUDE "bp.set_rss.f90"

!INCLUDE "efs.get.f90"



END MODULE es