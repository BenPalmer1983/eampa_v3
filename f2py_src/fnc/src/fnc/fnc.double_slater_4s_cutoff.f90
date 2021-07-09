!# Slater 4S

SUBROUTINE double_slater_4s_cutoff(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
IF(r .GT. p_fixed(1))THEN
  y = 0.0D0
ELSE
  y = (p_fixed(1) - r)**3 * ((p(1) * r**3 * exp(-1.0D0 * p(2) * r))**2 + (p(3) * r**3 * exp(-1.0D0 * p(4) * r)))**2
END IF
END SUBROUTINE double_slater_4s_cutoff





SUBROUTINE double_slater_4s_cutoff_v(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL double_slater_4s_cutoff(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE double_slater_4s_cutoff_v
