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
USE expspline, ONLY: exp_spline_fit_a, exp_spline_fit_b

IMPLICIT NONE

INCLUDE "fnc.globals.f90"
INCLUDE "fnc.interfaces.f90"

CONTAINS


INCLUDE "fnc.f.f90"
INCLUDE "fnc.fv.f90"
INCLUDE "fnc.fgrad.f90"
INCLUDE "fnc.str_compare.f90"
INCLUDE "fnc.fzbl.f90"
INCLUDE "fnc.interp.f90"
INCLUDE "fnc.step.f90"
INCLUDE "fnc.heaviside.f90"
INCLUDE "fnc.generic.f90"
INCLUDE "fnc.minverse.f90"


! AVAILABLE FUNCTIONS
INCLUDE "fnc.ackland_embedding.f90"
INCLUDE "fnc.ackland_embedding_grad.f90"
INCLUDE "fnc.ackland_mendelev_pair.f90"
INCLUDE "fnc.ackland_mendelev_pair_grad.f90"
INCLUDE "fnc.ackland_mendelev_pair_zbl.f90"
!INCLUDE "fnc.ackland_mendelev_pair_zbl_grad.f90"
INCLUDE "fnc.buckingham.f90"
INCLUDE "fnc.buckingham_grad.f90"
INCLUDE "fnc.cubic_knot_spline.f90"
INCLUDE "fnc.cubic_knot_spline_grad.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end_grad.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end_pair.f90"
INCLUDE "fnc.cubic_knot_spline_fixed_end_pair_grad.f90"
INCLUDE "fnc.cubic_knot_spline_2.f90"
INCLUDE "fnc.cubic_knot_spline_2_grad.f90"
INCLUDE "fnc.cubic_knot_spline_3.f90"
INCLUDE "fnc.cubic_knot_spline_3_grad.f90"
INCLUDE "fnc.cubic_knot_spline_4.f90"
INCLUDE "fnc.cubic_knot_spline_4_grad.f90"
INCLUDE "fnc.cubic_knot_spline_5.f90"
INCLUDE "fnc.cubic_knot_spline_5_grad.f90"
INCLUDE "fnc.cubic_spline.f90"
INCLUDE "fnc.cubic_spline_grad.f90"
INCLUDE "fnc.cubic_spline_zbl.f90"
INCLUDE "fnc.cubic_spline_zbl_grad.f90"
INCLUDE "fnc.cubic_spline_zbl_2.f90"
INCLUDE "fnc.cubic_spline_zbl_2_grad.f90"
INCLUDE "fnc.double_slater_4s.f90"
INCLUDE "fnc.double_slater_4s_grad.f90"
INCLUDE "fnc.double_slater_4s_cutoff.f90"
INCLUDE "fnc.double_slater_4s_cutoff_grad.f90"
INCLUDE "fnc.embedding_a.f90"
INCLUDE "fnc.embedding_a_grad.f90"
INCLUDE "fnc.embedding_b.f90"
INCLUDE "fnc.embedding_b_grad.f90"
INCLUDE "fnc.embedding_c.f90"
INCLUDE "fnc.embedding_c_grad.f90"
INCLUDE "fnc.embedding_d.f90"
INCLUDE "fnc.embedding_d_grad.f90"
INCLUDE "fnc.embedding_e.f90"
INCLUDE "fnc.embedding_e_grad.f90"
INCLUDE "fnc.embedding_f.f90"
INCLUDE "fnc.embedding_f_grad.f90"
INCLUDE "fnc.embedding_g.f90"
INCLUDE "fnc.embedding_g_grad.f90"
INCLUDE "fnc.embedding_h.f90"
INCLUDE "fnc.embedding_h_grad.f90"
INCLUDE "fnc.fs_embedding.f90"
INCLUDE "fnc.fs_embedding_grad.f90"
INCLUDE "fnc.lennard_jones.f90"
INCLUDE "fnc.lennard_jones_grad.f90"
INCLUDE "fnc.mendelev_embedding.f90"
INCLUDE "fnc.mendelev_embedding_grad.f90"
INCLUDE "fnc.morse.f90"
INCLUDE "fnc.morse_grad.f90"
INCLUDE "fnc.quad_embedding.f90"
INCLUDE "fnc.quad_embedding_grad.f90"
INCLUDE "fnc.quadratic_density.f90"
INCLUDE "fnc.quadratic_density_grad.f90"
INCLUDE "fnc.slater_4s.f90"
INCLUDE "fnc.slater_4s_grad.f90"
INCLUDE "fnc.slater_4s_cutoff.f90"
INCLUDE "fnc.slater_4s_cutoff_grad.f90"
INCLUDE "fnc.zbl.f90"
INCLUDE "fnc.zbl_grad.f90"
INCLUDE "fnc.zero.f90"
INCLUDE "fnc.zero_grad.f90"







! PAIR
INCLUDE "fnc.quartic_poly_rep.f90"
INCLUDE "fnc.pair_spline.f90"
INCLUDE "fnc.triple_embedding.f90"


! GENERAL FUNCTIONS ANY TIME
INCLUDE "fnc.quartic_poly.f90"
INCLUDE "fnc.quintic_spline.f90"

! OLD
INCLUDE "fnc.spline_ab.f90"
INCLUDE "fnc.spline_one_node.f90"
INCLUDE "fnc.spline_n_node.f90"


END MODULE fnc
