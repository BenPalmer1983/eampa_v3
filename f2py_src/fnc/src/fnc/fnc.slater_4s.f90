!# Slater 4S

SUBROUTINE slater_4s(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = (p(1) * r**3 * exp(-1.0D0 * p(2) * r))**2
END SUBROUTINE slater_4s


SUBROUTINE slater_4s_v(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL slater_4s(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE slater_4s_v
