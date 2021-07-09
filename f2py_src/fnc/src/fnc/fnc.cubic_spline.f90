! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
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
REAL(kind=DoubleReal) :: H, rk
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
DO n = 1, SIZE(p,1)
  rk = p_fixed(n)
  CALL heaviside((rk - r), H)
  y = y + H * p(n) * (r - rk)**3
END DO
END SUBROUTINE cubic_spline

! VECTOR SUBROUTINE
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
! Loop through all the values in r(:), calculate and store in y(:)
DO n = 1, SIZE(r,1)
  CALL cubic_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_v
