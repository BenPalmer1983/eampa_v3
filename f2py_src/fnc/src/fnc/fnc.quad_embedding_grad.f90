SUBROUTINE quad_embedding_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = A * sqrt(r) + B * r + C * r**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: h = 1.0D-7
REAL(kind=DoubleReal) :: a, b
!############################################################
!CALL fs_embedding(r-h, p, p_fixed, a)
!CALL fs_embedding(r+h, p, p_fixed, b)
!dydr = (b - a) / (2.0D0 * h)
dydr = 0.5 * p(2) * r**(-0.5) + 2 * p(3) * r + 4 * p(4) * r**3
END SUBROUTINE quad_embedding_grad