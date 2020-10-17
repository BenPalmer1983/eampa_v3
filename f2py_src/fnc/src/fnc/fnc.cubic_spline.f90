! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

SUBROUTINE cubic_spline(r, p, p_fixed, y)
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
  y = y + p(n) * (r - p_fixed(n))**3 * H
END DO
END SUBROUTINE cubic_spline

SUBROUTINE cubic_spline_v(r, p, p_fixed, y)
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
  CALL cubic_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_v
