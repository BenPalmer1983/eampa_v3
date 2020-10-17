!#############################################################
!#  HOW TO USE:
!#
!#
!#                                                        
!#
!#############################################################

MODULE fnc

USE kinds
USE sls, ONLY: sls_solve

IMPLICIT NONE

INCLUDE "fnc.globals.f90"
INCLUDE "fnc.interfaces.f90"

CONTAINS

! STEP FUNCTIONS
INCLUDE "fnc.step.f90"
INCLUDE "fnc.heaviside.f90"

! PAIR
INCLUDE "fnc.lennard_jones.f90"
INCLUDE "fnc.morse.f90"
INCLUDE "fnc.buckingham.f90"
INCLUDE "fnc.zbl.f90"
INCLUDE "fnc.quartic_poly_rep.f90"

! DENS
INCLUDE "fnc.quadratic_density.f90"
INCLUDE "fnc.slater_4s.f90"

!EMBED
INCLUDE "fnc.fs_embedding.f90"
INCLUDE "fnc.mendelev_embedding.f90"
INCLUDE "fnc.triple_embedding.f90"
INCLUDE "fnc.ackland_embedding.f90"


! GENERAL FUNCTIONS ANY TIME
INCLUDE "fnc.quartic_poly.f90"
INCLUDE "fnc.cubic_spline.f90"
INCLUDE "fnc.quintic_spline.f90"




! NODE TO NODE SPLINES
INCLUDE "fnc.spline_ab.f90"
INCLUDE "fnc.spline_one_node.f90"
INCLUDE "fnc.spline_n_node.f90"


END MODULE fnc
