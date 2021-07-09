

SUBROUTINE ackland_embedding_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: h
REAL(kind=DoubleReal) :: a, b
!############################################################
h = 1.0D-7
CALL ackland_embedding(r-h, p, p_fixed, a)
CALL ackland_embedding(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE ackland_embedding_grad


