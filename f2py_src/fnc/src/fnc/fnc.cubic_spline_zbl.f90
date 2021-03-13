! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
SUBROUTINE cubic_spline_zbl(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: Q1, Q2, rc1, rc2, rk, rs, x, e, ep
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb, H
INTEGER(kind=StandardInteger) :: n_p, n_pf
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
!############################################################
Q1 = p_fixed(1)
Q2 = p_fixed(2)
rc1 = p_fixed(3)


IF(r .EQ. 0.0D0)THEN
  y = 1.0D4
ELSE
  CALL heaviside(rc1 - r, H)  
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (r / rs)
  e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
  y = H * (rc1 - r) * ((Q1 * Q2) / r) * e
  DO n_p = 1, SIZE(p,1)
    n_pf = n_p + 3
    rk = p_fixed(n_pf)
    CALL heaviside(r - rk, H)
    y = y + p(n_p) * H * (r - rk)**3
  END DO
END IF


END SUBROUTINE cubic_spline_zbl



! VECTOR SUBROUTINE
SUBROUTINE cubic_spline_zbl_v(r, p, p_fixed, y)
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
  CALL cubic_spline_zbl(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_zbl_v


