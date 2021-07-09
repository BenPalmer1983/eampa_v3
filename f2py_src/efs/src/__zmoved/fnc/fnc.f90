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



INCLUDE "fnc.fzbl.f90"
INCLUDE "fnc.interp.f90"

! STEP FUNCTIONS
INCLUDE "fnc.step.f90"
INCLUDE "fnc.heaviside.f90"

! PAIR
INCLUDE "fnc.lennard_jones.f90"
INCLUDE "fnc.morse.f90"
INCLUDE "fnc.buckingham.f90"
INCLUDE "fnc.zbl.f90"
INCLUDE "fnc.quartic_poly_rep.f90"
INCLUDE "fnc.ackland_mendelev_pair.f90"
INCLUDE "fnc.cubic_spline_zbl.f90"
INCLUDE "fnc.pair_spline.f90"



! DENS
INCLUDE "fnc.quadratic_density.f90"
INCLUDE "fnc.slater_4s.f90"
INCLUDE "fnc.double_slater_4s_cutoff.f90"


!EMBED
INCLUDE "fnc.fs_embedding.f90"
INCLUDE "fnc.mendelev_embedding.f90"
INCLUDE "fnc.triple_embedding.f90"
INCLUDE "fnc.quad_embedding.f90"
INCLUDE "fnc.ackland_embedding.f90"


! GENERAL FUNCTIONS ANY TIME
INCLUDE "fnc.quartic_poly.f90"
INCLUDE "fnc.cubic_spline.f90"
INCLUDE "fnc.quintic_spline.f90"




! NODE TO NODE SPLINES
INCLUDE "fnc.cubic_knot_spline.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end_pair.f90"



! OLD
INCLUDE "fnc.spline_ab.f90"
INCLUDE "fnc.spline_one_node.f90"
INCLUDE "fnc.spline_n_node.f90"


END MODULE fnc
