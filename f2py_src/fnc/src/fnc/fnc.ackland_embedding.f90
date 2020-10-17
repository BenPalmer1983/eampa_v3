!# See Olsson, Wallenius 2005

SUBROUTINE ackland_embedding(r, p, p_fixed, y)
!############################################################
! f(x) = A * sqrt(r) + B * r + C * r**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = p(1) * sqrt(r) + p(2) * r**2 + p(3) * r**4
END SUBROUTINE ackland_embedding


SUBROUTINE ackland_embedding_v(r, p, p_fixed, y)
!############################################################
! f(x) = A * sqrt(r) + B * r + C * r**2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL ackland_embedding(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE ackland_embedding_v
