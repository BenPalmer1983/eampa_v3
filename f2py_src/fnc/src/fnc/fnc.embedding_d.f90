! Development of an interatomic potential for phosphorus impurities in α-iron
! 2004

SUBROUTINE embedding_d(r, p, p_fixed, y)
!############################################################
! f(x) = -1.0D0 * sqrt(r) + p(1) * r**2 + p(2) * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
IF(r .LT. 0.0D0)THEN
  y = 0.0D0
ELSE
  y = -1.0D0 * sqrt(r) + p(1) * r**2 + p(2) * r**4
END IF
END SUBROUTINE embedding_d



!Development of an interatomic potential for phosphorus impurities in α-iron


SUBROUTINE embedding_d_v(r, p, p_fixed, y)
!############################################################
! f(x) = -sqrt(rho) + A*rho**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL embedding_d(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE embedding_d_v