


SUBROUTINE spline_ab(node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: node_a(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: node_b(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: coeffs(1:6)
!############################################################
INTEGER(kind=StandardInteger) :: m_size = 0
INTEGER(kind=StandardInteger) :: m_half_size = 0
INTEGER(kind=StandardInteger) :: d_loop = 0
REAL(kind=DoubleReal) :: x(1:6,1:6)
REAL(kind=DoubleReal) :: y(1:6)
REAL(kind=DoubleReal) :: x_coeffs(1:6)
REAL(kind=DoubleReal) :: x_exponent(1:6)
INTEGER(kind=StandardInteger) :: n, row, col, a_row
!############################################################
! Sizes
m_size = 6
m_half_size = 3
! x exponents
DO n = 1, m_size
  x_exponent(n) = 1.0D0 * (n - 1)
END DO
! Coeffs
x_coeffs = 1
! Make y matrix
row = 1
DO n=2,4
  y(row) = node_a(n)
  row = row + 1
  y(row) = node_b(n)
  row = row + 1
END DO
! Make x matrix
row = 0
DO a_row = 1, m_half_size
  ! NODE A
  row = row + 1
  DO col = 1, a_row -1
    x(row, col) = 0.0D0
  END DO
  DO col = a_row, m_size
    x(row, col) = x_coeffs(col) * node_a(1)**x_exponent(col)
  END DO
  ! NODE B
  row = row + 1
  DO col = 1, a_row -1
    x(row, col) = 0.0D0
  END DO
  DO col = a_row, m_size
    x(row, col) = x_coeffs(col) * node_b(1)**x_exponent(col)
  END DO
  ! UPDATE COEFFS AND EXPONENT
  DO col = 1, m_size
    x_coeffs(col) = x_coeffs(col) * x_exponent(col)
  END DO
  DO col = 1, m_size
    IF(x_exponent(col) .GT. 0.0D0)THEN
      x_exponent(col) = x_exponent(col) - 1.0D0
    END IF
  END DO
END DO
! SOLVE
CALL sls_solve(x, y, coeffs)
!############################################################
END SUBROUTINE spline_ab
