

SUBROUTINE poly_row(coeffs, x, row)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: coeffs(:)
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(INOUT) :: row(:)
!############################################################
INTEGER(kind=StandardInteger) :: n, i , j
REAL(kind=DoubleReal) :: c(1:SIZE(coeffs))
!############################################################
c(:) = coeffs(:)
row(:) = 0.0D0
row(1) = x
DO i=2,SIZE(row)
  DO j=1,SIZE(c,1)
    row(i) = row(i) + c(j) * x**(j-i+1)
  END DO
  DO j=1, SIZE(c,1)
    c(j) = (j-i+1)*c(j)
  END DO
END DO
END SUBROUTINE poly_row



SUBROUTINE poly_exp_row(coeffs, x, row)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: coeffs(:)
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(INOUT) :: row(:)
!############################################################
INTEGER(kind=StandardInteger) :: n, i , j
REAL(kind=DoubleReal) :: c(1:SIZE(coeffs))
!############################################################
c(:) = coeffs(:)
row(:) = 0.0D0
row(1) = x
DO i=2,SIZE(row)
  DO j=1,SIZE(c,1)
    row(i) = row(i) + c(j) * x**(j-i+1)
  END DO
  row(i) = EXP(row(i))
  DO j=1, SIZE(c,1)
    c(j) = (j-i+1)*c(j)
  END DO
END DO
END SUBROUTINE poly_exp_row






