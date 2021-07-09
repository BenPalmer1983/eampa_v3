SUBROUTINE zero_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = A * sqrt(r) + B * r**2 + C * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
!############################################################
dydr = 0.0D0
END SUBROUTINE zero_grad