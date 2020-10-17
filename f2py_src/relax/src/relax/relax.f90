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

MODULE relax

USE kinds
USE matrix, ONLY: inv
USE polyfit, ONLY: fit
USE interp, ONLY: fill, fill_zoor, interp4dydxn
USE spline, ONLY: spline_array, spline_ab_poly, spline_ab_exp3, spline_ab_exp5
USE sls, ONLY: solve

IMPLICIT NONE

INCLUDE "relax.globals.f90"
INCLUDE "relax.interfaces.f90"

CONTAINS

INCLUDE "relax.init.f90"
INCLUDE "relax.math.f90"

INCLUDE "relax.potential.f90"
INCLUDE "relax.add_config.f90"
INCLUDE "relax.refresh_coords.f90"
INCLUDE "relax.make_ghost.f90"
INCLUDE "relax.ghost_update.f90"
INCLUDE "relax.nl.f90"
INCLUDE "relax.nl_update.f90"
INCLUDE "relax.zero.f90"
INCLUDE "relax.force_calc.f90"
INCLUDE "relax.run.f90"
INCLUDE "relax.run_md.f90"



!INCLUDE "es.add_calculation.f90"


!INCLUDE "es.nl.f90"
!INCLUDE "es.energy.f90"
!INCLUDE "es.force.f90"


!INCLUDE "es.es.f90"


!INCLUDE "bp.add_config.f90"
!INCLUDE "bp.add_calculation.f90"
!INCLUDE "bp.add_known.f90"
!INCLUDE "bp.nl.f90"

!INCLUDE "bp.energy.f90"

!INCLUDE "bp.bp.f90"

!INCLUDE "bp.set_rss.f90"

!INCLUDE "efs.get.f90"



END MODULE relax