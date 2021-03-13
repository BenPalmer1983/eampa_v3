

!############################################################
! Fix the end nodes with p_fix (x, y, y', y'' for start and end)


SUBROUTINE spline_n_node(r, p_var, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p_var(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: n, c
REAL(kind=DoubleReal) :: coeffs(1:6)
!############################################################
n = 1
DO WHILE(n .LT. (SIZE(p_var, 1) / 4))
  IF(r .GE. p_var(4 * (n - 1) + 1) .AND. r .LE. p_var(4 * n + 1))THEN
    CALL spline_ab(p_var(4 * (n - 1) + 1:4 * (n - 1) + 4), p_var(4 * n + 1:4 * n + 4), coeffs(1:6))
    n = SIZE(p_var, 1)
  END IF
  n = n + 1
END DO
y = 0.0D0
DO c = 1, 6
  y = y + coeffs(c) * r**(c-1)
END DO
END SUBROUTINE spline_n_node




SUBROUTINE spline_n_node_v(r, p_var, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_var(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL spline_n_node(r(n), p_var(:), y(n))
END DO
END SUBROUTINE spline_n_node_v

