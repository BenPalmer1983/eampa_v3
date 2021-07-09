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

MODULE ce

USE kinds
USE potential, ONLY: f, search_f, search_grad, f_groups, fgroup_max
USE matrix, ONLY: inv
USE polyfit, ONLY: fit
USE sls, ONLY: solve

IMPLICIT NONE

INCLUDE "ce.globals.f90"
INCLUDE "ce.interfaces.f90"

CONTAINS

INCLUDE "ce.init.f90"
INCLUDE "ce.math.f90"



!INCLUDE "bp.add_config.f90"
!INCLUDE "bp.add_calculation.f90"
!INCLUDE "bp.add_known.f90"
!INCLUDE "bp.nl.f90"

!INCLUDE "bp.energy.f90"

!INCLUDE "bp.calculate_bp.f90"

!INCLUDE "bp.set_rss.f90"


!INCLUDE "bp.potential.f90"
!INCLUDE "efs.get.f90"



END MODULE ce