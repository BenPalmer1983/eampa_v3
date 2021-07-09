! VECTOR SUBROUTINE
SUBROUTINE fgrad(f_name, r, p, p_fixed, dydr)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "ackland_mendelev_pair"))THEN
  CALL ackland_mendelev_pair_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline"))THEN
  CALL cubic_knot_spline_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end"))THEN
  CALL cubic_knot_spline_fixed_end_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end_pair"))THEN
  CALL cubic_knot_spline_fixed_end_pair_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_2"))THEN
  CALL cubic_knot_spline_2_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_3"))THEN
  CALL cubic_knot_spline_3_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_4"))THEN
  CALL cubic_knot_spline_4_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_knot_spline_5"))THEN
  CALL cubic_knot_spline_5_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_spline"))THEN
  CALL cubic_spline_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_spline_zbl"))THEN
  CALL cubic_spline_zbl_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "cubic_spline_zbl_2"))THEN
  CALL cubic_spline_zbl_2_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "double_slater_4s"))THEN
  CALL double_slater_4s_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "double_slater_4s_cutoff"))THEN
  CALL double_slater_4s_cutoff_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_a"))THEN
  CALL embedding_a_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_b"))THEN
  CALL embedding_b_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_c"))THEN
  CALL embedding_c_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_d"))THEN
  CALL embedding_d_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_e"))THEN
  CALL embedding_e(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_f"))THEN
  CALL embedding_f(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_g"))THEN
  CALL embedding_g(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "embedding_h"))THEN
  CALL embedding_h(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "quad_embedding"))THEN
  CALL quad_embedding_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "quadratic_density"))THEN
  CALL quadratic_density_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "slater_4s_cutoff"))THEN
  CALL slater_4s_cutoff_grad(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl_grad(r, p, p_fixed, dydr)
END IF

END SUBROUTINE fgrad
