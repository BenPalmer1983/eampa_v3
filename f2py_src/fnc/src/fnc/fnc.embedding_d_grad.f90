SUBROUTINE embedding_d_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: h = 1.0D-7
REAL(kind=DoubleReal) :: a, b
!############################################################
CALL embedding_d(r-h, p, p_fixed, a)
CALL embedding_d(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE embedding_d_grad