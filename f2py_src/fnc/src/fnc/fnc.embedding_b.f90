! Development of new interatomic potentials appropriate for crystalline and liquid iron
! 2003

SUBROUTINE embedding_b(r, p, p_fixed, y)
!############################################################
! f(x) = -1.0D0 * sqrt(r) + p(1) * r**2
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
  y = -1.0D0 * sqrt(r) + p(1) * r**2
END IF
END SUBROUTINE embedding_b



SUBROUTINE embedding_b_v(r, p, p_fixed, y)
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
  CALL embedding_b(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE embedding_b_v