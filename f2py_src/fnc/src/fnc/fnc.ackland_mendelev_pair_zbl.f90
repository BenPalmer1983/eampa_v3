! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
SUBROUTINE ackland_mendelev_pair_zbl(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: Q1, Q2, rc1, rc2, rc3, rk, rs, x, e, ep, psize, b0, b1, b2, b3
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb, H1, H2
INTEGER(kind=StandardInteger) :: n_p, n_pf
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
REAL(kind=DoubleReal) :: pf(1:SIZE(p,1))
!############################################################

psize = SIZE(p,1)
IF(SIZE(p_fixed,1) .NE. SIZE(p,1) + 8)THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = 0.0D0
ELSE 

pf(1:psize) = p_fixed(1:psize)
Q1 = p_fixed(psize+1)
Q2 = p_fixed(psize+2)
rc1 = p_fixed(psize+3)
rc2 = p_fixed(psize+4)
b0 = p_fixed(psize+5) 
b1 = p_fixed(psize+6) 
b2 = p_fixed(psize+7) 
b3 = p_fixed(psize+8) 
rc3 = p_fixed(psize)   ! the last cutoff

IF(r .EQ. 0.0D0)THEN
  y = 1.0D4
ELSE IF(r .LE. 0.0D0 .OR. r .GT. rc3)THEN
  y = 0.0D0
ELSE
  y = 0.0D0
  
  ! ZBL
  CALL heaviside(rc1 - r, H1)
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (r / rs)
  e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
  y = H1 * ((Q1 * Q2) / r) * e

  ! SPLINE
  CALL heaviside(rc2 - r, H2)
	CALL heaviside(r - rc1, H1)
  y = y + H1 * H2 * exp(b0 + b1 * r + b2 * r**2 + b3 * r**3)
  
  ! Cubic spline
  DO n_p = 1, SIZE(p,1)
  	rk = p_fixed(n_pf)
    CALL heaviside(rk - r, H2)
	  CALL heaviside(r - rc2, H1)
    y = y + H1 * H2 * p(n_p) * (r - rk)**3
  END DO

END IF


END IF

END SUBROUTINE ackland_mendelev_pair_zbl







! VECTOR SUBROUTINE
SUBROUTINE ackland_mendelev_pair_zbl_v(r, p, p_fixed, y)
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
  CALL ackland_mendelev_pair_zbl(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE ackland_mendelev_pair_zbl_v


