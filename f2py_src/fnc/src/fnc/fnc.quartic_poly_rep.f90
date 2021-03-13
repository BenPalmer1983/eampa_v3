! V-Ti-Cr Fu, Li, Johansson, Zhao

SUBROUTINE quartic_poly_rep(r, p, p_fixed, y)
!############################################################
! Quartic Poly (r - rc)^2(c0 + c1 r + c2 r^2) + D exp (-(r/q))
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:5)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)  ! rcut
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
IF(r .LE. p_fixed(1))THEN
  y = (r - p_fixed(1))**2 * (p(1) + p(2) * r + p(3) * r**2) + p(4) * exp(-1.0D0 * (r / p(5)))
ELSE
  y = 0.0D0
END IF
END SUBROUTINE quartic_poly_rep

SUBROUTINE quartic_poly_rep_v(r, p, p_fixed, y)
!############################################################
! Quartic Poly (r - rc)^2(c0 + c1 r + c2 r^2)
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:5)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL quartic_poly_rep(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE quartic_poly_rep_v