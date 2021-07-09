
SUBROUTINE embedding_e(r, p, p_fixed, y)
!############################################################
! f(x) = p(1) * sqrt(r) + p(2) * r**2 + p(3) * r**4
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
  y = p(1) * sqrt(r) + p(2) * r**2 + p(3) * r**4
END IF
END SUBROUTINE embedding_e



SUBROUTINE embedding_e_v(r, p, p_fixed, y)
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
  CALL embedding_e(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE embedding_e_v