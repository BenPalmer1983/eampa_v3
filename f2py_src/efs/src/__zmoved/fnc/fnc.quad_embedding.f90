!# See Setyawan, Gao, Kurtz W-Re 2018

SUBROUTINE quad_embedding(r, p, p_fixed, y)
!############################################################
! f(x) = A * sqrt(r) + B * r + C * r**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = p(1) + p(2) * sqrt(r) + p(3) * r**2 + p(4) * r**4
END SUBROUTINE quad_embedding

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
!############################################################
dydr = 0.5 * p(2) * r**(-0.5) + 2 * p(3) * r + 4 * p(4) * r**3
END SUBROUTINE quad_embedding_grad



SUBROUTINE quad_embedding_v(r, p, p_fixed, y)
!############################################################
! f(x) = A * sqrt(r) + B * r + C * r**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL quad_embedding(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE quad_embedding_v
