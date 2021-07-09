SUBROUTINE pivot(A, A_out)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: A(:,:)
REAL(kind=DoubleReal), INTENT(OUT) :: A_out(1:SIZE(A,1),1:SIZE(A,2))
!############################################################
REAL(kind=DoubleReal) :: PivotMatrix(1:SIZE(A,1),1:SIZE(A,2))
INTEGER(kind=StandardInteger) :: n, j
INTEGER(kind=StandardInteger) :: max_row
REAL(kind=DoubleReal) :: max_row_val
REAL(kind=DoubleReal) :: temp
!############################################################
!# Get identity
PivotMatrix = 0.0D0
DO j = 1, SIZE(A,1)
  PivotMatrix(j,j) = 1.0D0
END DO
DO j = 1, SIZE(A,1)
  max_row = j
  max_row_val = 0.0D0
  DO n = j, SIZE(A,2)
    IF (ABS(A(n,j)) .GT. max_row_val) THEN
      max_row_val = ABS(A(n,j))
      max_row = n
    END IF
  END DO
  IF (j .NE. max_row) THEN
    !# Swap row
    DO n = 1, SIZE(A,2)
      temp = PivotMatrix(j,n)
      PivotMatrix(j,n) = PivotMatrix(max_row,n)
      PivotMatrix(max_row,n) = temp
    END DO
  END IF
END DO
A_out = MATMUL(PivotMatrix, A)
!############################################################
END SUBROUTINE pivot