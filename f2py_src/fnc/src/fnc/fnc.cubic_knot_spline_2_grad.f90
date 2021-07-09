SUBROUTINE cubic_knot_spline_2_grad(r, p, p_fixed, dydr)
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
CALL cubic_knot_spline_2(r-h, p, p_fixed, a)
CALL cubic_knot_spline_2(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE cubic_knot_spline_2_grad