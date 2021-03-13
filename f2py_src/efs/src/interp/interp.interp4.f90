SUBROUTINE interp4(xi, x, y, yi)
! Identity for square matrix
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: y(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: yi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j
REAL(kind=DoubleReal) :: li
!############################################################
yi = 0.0D0
IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
  DO i = 1, 4
    li = 1.0D0
    DO j = 1, 4
      IF(i .NE. j) THEN
        li = li * (xi - x(j)) / (x(i) - x(j))
      END IF
    END DO
    yi = yi + li * y(i)
  END DO
END IF
!############################################################
END SUBROUTINE interp4