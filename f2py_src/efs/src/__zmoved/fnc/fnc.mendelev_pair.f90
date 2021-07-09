! Mendelev Pair
! Ackland, Mendelev, Srolovitz, Han, Barashev 2004
! MENDELEV_PAIR   ZBL + CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

SUBROUTINE mendelev_pair(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: qi, qj, za, zb
REAL(kind=DoubleReal) :: zbl
REAL(kind=DoubleReal) :: rs, rr
REAL(kind=DoubleReal) :: a0, b0, b1, b2, b3
REAL(kind=DoubleReal) :: H
INTEGER(kind=StandardInteger) :: n
!############################################################
qi = p_fixed(1)
qj = p_fixed(2)
za = p_fixed(3)
zb = p_fixed(4)
a0 = p(1)
b0 = p(2)
b1 = p(3)
b2 = p(4)
b3 = p(5)

IF(r .EQ. 0.0D0)THEN
  y = 1.0D99
ELSE IF(r .GT. 0.0D0 .AND. r .LE. za)THEN
  rs = 0.46837766 / (qi**(2.0D0/3.0D0) + qj**(2.0D0/3.0D0))
  rr = r / rs
  zbl = 0.1818D0 * exp(-3.2D0 * rr) + 0.5099D0 * exp(-0.9423D0 * rr) + 0.2802 * exp(-0.4029 * rr) + 0.02817 * exp(-0.2016 * rr)
  y = a0 * ((qi * qi) / r) * zbl
ELSE IF(r .GT. p_fixed(3) .AND. r .LE. zb)THEN
  y = exp(b0 + b1 * r + b2 * r**2 + b3 * r**3)
ELSE IF(r .GT. zb)THEN
  y = 0.0D0
  DO n = 1, SIZE(p,1)-6
    IF(r .LT. p_fixed(4+n))THEN
      y = y + p(n) * (r - p_fixed(4 + n))**3      
    END IF
  END DO
END IF


!  IF((p_fixed(n+6) - r) < 0.0D0)THEN
!    H = 0.0D0
!  ELSE
!    H = 1.0D0
!  END IF
!  y = y + p(n) * (r - p_fixed(n))**3 * H
!END DO
END SUBROUTINE mendelev_pair

SUBROUTINE mendelev_pair_v(r, p, p_fixed, y)
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
  CALL mendelev_pair(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE mendelev_pair_v