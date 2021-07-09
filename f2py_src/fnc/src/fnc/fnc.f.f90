! VECTOR SUBROUTINE
SUBROUTINE f(f_name, r, p, p_fixed, y)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_mendelev_pair"))THEN
  CALL ackland_mendelev_pair(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_mendelev_pair_zbl"))THEN
  CALL ackland_mendelev_pair_zbl(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline"))THEN
  CALL cubic_knot_spline(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end"))THEN
  CALL cubic_knot_spline_fixed_end(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_fixed_end_pair"))THEN
  CALL cubic_knot_spline_fixed_end_pair(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_2"))THEN
  CALL cubic_knot_spline_2(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_3"))THEN
  CALL cubic_knot_spline_3(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_4"))THEN
  CALL cubic_knot_spline_4(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_knot_spline_5"))THEN
  CALL cubic_knot_spline_5(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline"))THEN
  CALL cubic_spline(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline_zbl"))THEN
  CALL cubic_spline_zbl(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "cubic_spline_zbl_2"))THEN
  CALL cubic_spline_zbl_2(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "double_slater_4s"))THEN
  CALL double_slater_4s(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "double_slater_4s_cutoff"))THEN
  CALL double_slater_4s_cutoff(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_a"))THEN
  CALL embedding_a(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_b"))THEN
  CALL embedding_b(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_c"))THEN
  CALL embedding_c(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_d"))THEN
  CALL embedding_d(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_e"))THEN
  CALL embedding_e(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_f"))THEN
  CALL embedding_f(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_g"))THEN
  CALL embedding_g(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "embedding_h"))THEN
  CALL embedding_h(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "quad_embedding"))THEN
  CALL quad_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "quadratic_density"))THEN
  CALL quadratic_density(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s_cutoff"))THEN
  CALL slater_4s_cutoff(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl(r, p, p_fixed, y)
END IF

END SUBROUTINE f
