! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
SUBROUTINE ackland_mendelev_pair(r, p, p_fixed, y)
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
rc2 = p_fixed(4)


IF(r .EQ. 0.0D0)THEN
  y = 1.0D4
ELSE IF(r .LE. rc1)THEN
  CALL heaviside(rc1 - r, H)
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (r / rs)
  e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
  y = H * ((Q1 * Q2) / r) * e
ELSE IF(r .GT. rc1 .AND. r .LT. rc2)THEN
  ! x points
  xa = rc1
  xb = rc2
  
  ! ya y'a
  CALL heaviside(xa - r, H)
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (xa / rs)
  e = 0.1818D0 * exp(-3.2D0 * x) + 0.5099D0 * exp(-0.9423D0 * x) + 0.2802D0 * exp(-0.4029D0 * x) + 0.02817D0 * exp(-0.2016D0 * x)
  ep =      (-3.2D0 / rs) * 0.1818D0  * exp(-3.2D0 * x)         + (-0.9423D0 / rs) * 0.5099D0 * exp(-0.9423D0 * x)
  ep = ep + (-0.4029D0 / rs) * 0.2802D0  * exp(-0.4029D0 * x)   + (-0.2016D0 / rs) * 0.02817D0  * exp(-0.2016D0 * x)
  ya = H * ((Q1 * Q2) / xa) * e
  ypa = -1.0D0 * H * ((Q1 * Q2)/ xa**2) * e + H * ((Q1 * Q2) / xa) * ep
  
  yb = 0.0D0
  DO n_p = 1, SIZE(p,1)
    rk = p_fixed(n_pf + 4)
    CALL heaviside(xb - rk, H)
    yb = yb + p(n_p) * H * (xb - rk)**3
  END DO
  
  ypb = 0.0D0
  DO n_p = 1, SIZE(p,1)
    n_pf = n_p + 4
    rk = p_fixed(n_pf)
    CALL heaviside(xb - rk, H)
    ypb = ypb + 3 * p(n_p) * H * (xb - rk)**2
  END DO
 
  
  xmat(1,1) = 1.0
  xmat(1,2) = xa
  xmat(1,3) = xa**2
  xmat(1,4) = xa**3
  xmat(2,1) = 1.0
  xmat(2,2) = xb
  xmat(2,3) = xb**2
  xmat(2,4) = xb**3
  xmat(3,1) = 0.0
  xmat(3,2) = 1.0D0
  xmat(3,3) = 2.0D0 * xa
  xmat(3,4) = 3.0D0 * xa**2
  xmat(4,1) = 0.0
  xmat(4,2) = 1.0D0
  xmat(4,3) = 2.0D0 * xb
  xmat(4,4) = 3.0D0 * xb**2
  
  ymat(1) = ya
  ymat(2) = ypa
  ymat(3) = yb
  ymat(4) = ypb
  
  CALL sls_solve(xmat, ymat, c)
  y = c(1) + c(2) * r + c(3) * r**2 + c(4) * r**3

ELSE IF(r .GE. rc2)THEN
  y = 0.0D0
  DO n_p = 1, SIZE(p,1)
    n_pf = n_p + 4
    rk = p_fixed(n_pf)
    CALL heaviside(r - rk, H)
    y = y + p(n_p) * H * (r - rk)**3
  END DO
END IF


END SUBROUTINE ackland_mendelev_pair



! VECTOR SUBROUTINE
SUBROUTINE ackland_mendelev_pair_v(r, p, p_fixed, y)
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
  CALL ackland_mendelev_pair(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE ackland_mendelev_pair_v


