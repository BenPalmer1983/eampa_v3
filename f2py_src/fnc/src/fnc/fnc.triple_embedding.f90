!# See Setyawan, Gao, Kurtz W-Re 2018

SUBROUTINE triple_embedding(r, p, p_fixed, y)
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
y = p(1) * sqrt(r) + p(2) * r + p(3) * r**2
END SUBROUTINE triple_embedding


SUBROUTINE triple_embedding_v(r, p, p_fixed, y)
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
  CALL triple_embedding(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE triple_embedding_v
