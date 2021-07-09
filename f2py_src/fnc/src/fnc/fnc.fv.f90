! VECTOR SUBROUTINE
SUBROUTINE fv(f_name, r, p, p_fixed, y)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_mendelev_pair"))THEN
  CALL ackland_mendelev_pair_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline"))THEN
  CALL cubic_knot_spline_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end"))THEN
  CALL cubic_knot_spline_fixed_end_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end_pair"))THEN
  CALL cubic_knot_spline_fixed_end_pair_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_2"))THEN
  CALL cubic_knot_spline_2_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_3"))THEN
  CALL cubic_knot_spline_3_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_4"))THEN
  CALL cubic_knot_spline_4_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_5"))THEN
  CALL cubic_knot_spline_5_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline"))THEN
  CALL cubic_spline_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline_zbl"))THEN
  CALL cubic_spline_zbl_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline_zbl_2"))THEN
  CALL cubic_spline_zbl_2_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "double_slater_4s"))THEN
  CALL double_slater_4s_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "double_slater_4s_cutoff"))THEN
  CALL double_slater_4s_cutoff_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_a"))THEN
  CALL embedding_a_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_b"))THEN
  CALL embedding_b_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_c"))THEN
  CALL embedding_c_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_d"))THEN
  CALL embedding_d_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_e"))THEN
  CALL embedding_e_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_f"))THEN
  CALL embedding_f_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_g"))THEN
  CALL embedding_g_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_h"))THEN
  CALL embedding_h_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "quad_embedding"))THEN
  CALL quad_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "quadratic_density"))THEN
  CALL quadratic_density_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s_cutoff"))THEN
  CALL slater_4s_cutoff_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl_v(r, p, p_fixed, y)
END IF

END SUBROUTINE fv
