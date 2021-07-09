! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

SUBROUTINE cubic_spline_cutoff(r, p, p_fixed, y)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: H
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
DO n = 1, SIZE(p,1)
  IF((p_fixed(n) - r) < 0.0D0)THEN
    H = 0.0D0
  ELSE
    H = 1.0D0
  END IF
  y = y + p(n) * (p_fixed(n) - r)**3 * H
END DO
IF(SIZE(p,1) + 1 .EQ. SIZE(p_fixed,1))THEN
  IF(r .GE. p_fixed(SIZE(p_fixed,1)))THEN
    y = 0.0D0
  ELSE
    y = y * (p_fixed(SIZE(p_fixed,1)) - r)**3
  END IF
END IF
END SUBROUTINE cubic_spline_cutoff

SUBROUTINE cubic_spline_cutoff_v(r, p, p_fixed, y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL cubic_spline_cutoff(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_cutoff_v
