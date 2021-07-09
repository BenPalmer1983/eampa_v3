

SUBROUTINE Modulus(x, divisor, y)
!###########################################################
Real(kind=DoubleReal), INTENT(IN) :: x
Real(kind=DoubleReal), INTENT(IN) :: divisor
Real(kind=DoubleReal), INTENT(OUT) :: y
Real(kind=DoubleReal) :: factor
!###########################################################
If(x.lt.0.0D0)Then
  factor = ceiling(abs(x/(1.0D0*divisor)))
  y = x+factor*divisor
Else
  factor = floor(x/(1.0D0*divisor))
  y = x-factor*divisor
End If
END SUBROUTINE Modulus








FUNCTION CrossProduct(a, b) RESULT (c)
! Declare variables
REAL(kind=DoubleReal), DIMENSION(1:3) :: a, b, c
! Calculate cross product
c(1) = a(2)*b(3)-a(3)*b(2)
c(2) = a(3)*b(1)-a(1)*b(3)
c(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION CrossProduct




FUNCTION DotProduct(a, b) RESULT (c)
! Declare variables
REAL(kind=DoubleReal), DIMENSION(1:3) :: a, b
REAL(kind=DoubleReal) :: c
! Calculate cross product
c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
END FUNCTION DotProduct




FUNCTION TripleProduct(a, b, c) RESULT (v)
! Calculates (scalar) triple product of three vectors (resulting in volume of 3 vectors)
REAL(kind=DoubleReal), Dimension(1:3) :: a, b, c
REAL(kind=DoubleReal) :: v
! Calculate cross product
v = DotProduct(a,CrossProduct(b,c))
END FUNCTION TripleProduct




SUBROUTINE poly_calc(coeffs, x, fx)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: coeffs(:)
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(INOUT) :: fx
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
fx = 0.0D0
DO n=1,SIZE(coeffs,1)
  fx = fx + coeffs(n) * x**(n-1)
END DO
END SUBROUTINE poly_calc




SUBROUTINE exp_poly_calc(coeffs, x, fx)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: coeffs(:)
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(INOUT) :: fx
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
fx = 0.0D0
DO n=1,SIZE(coeffs,1)
  fx = fx + coeffs(n) * x**(n-1)
END DO
fx = EXP(fx)
END SUBROUTINE exp_poly_calc






