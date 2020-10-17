! SPLINE BETWEEN TWO NODES


! EXAMPLE 5th poly
!
! ALLOCATE(node_a(1:4))
! ALLOCATE(node_b(1:4))
! node_a(1) = 3
! node_a(2) = 3
! node_a(3) = 4
! node_a(4) = 1
! node_b(1) = 9
! node_b(2) = 6
! node_b(3) = -1
! node_b(4) = 1
! CALL spline_ab_poly(node_a, node_b, coeffs)
!



SUBROUTINE spline_ab(spline_type, node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: spline_type
REAL(kind=DoubleReal), INTENT(IN) :: node_a(:)
REAL(kind=DoubleReal), INTENT(IN) :: node_b(:)
REAL(kind=DoubleReal), INTENT(OUT) :: coeffs(1:SIZE(node_a) + SIZE(node_b) - 2)
!############################################################

coeffs = 0.0D0
IF(spline_type .EQ. 1 .AND. SIZE(node_a) .EQ. 3 .AND. SIZE(node_b) .EQ. 3)THEN
  CALL spline_ab_poly(node_a, node_b, coeffs)
END IF




END SUBROUTINE spline_ab




SUBROUTINE spline_ab_poly(node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: node_a
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: node_b
REAL(kind=DoubleReal), INTENT(OUT) :: coeffs(1:SIZE(node_a) + SIZE(node_b) - 2)
!############################################################
INTEGER(kind=StandardInteger) :: m_size = 0
INTEGER(kind=StandardInteger) :: m_half_size = 0
INTEGER(kind=StandardInteger) :: d_loop = 0
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:,:) :: x
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:) :: y
INTEGER(kind=StandardInteger) :: n, row, col, a_row
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:) :: x_coeffs
!############################################################
m_size = SIZE(node_a) + SIZE(node_b) - 2
m_half_size = m_size / 2
d_loop = SIZE(node_a) - 1

ALLOCATE(x(1:m_size,1:m_size))
ALLOCATE(x_coeffs(1:m_size))
ALLOCATE(y(1:m_size))

! Make y matrix
row = 1
DO n=2,SIZE(node_a)
  y(row) = node_a(n)
  row = row + 1
  y(row) = node_b(n)
  row = row + 1
END DO
! Make x matrix
row = 1
x_coeffs = 1
DO a_row = 1, d_loop
  ! NODE A
  DO col = 1, m_size 
    x(row, col) = x_coeffs(col) * node_a(1)**(col - a_row)
  END DO
  row = row + 1
  ! NODE B
  DO col = 1, m_size 
    x(row, col) = x_coeffs(col) * node_b(1)**(col - a_row)
  END DO
  row = row + 1
  ! UPDATE x_coeffs
  DO col = 1, m_size 
    x_coeffs(col) = x_coeffs(col) * (col - a_row)
  END DO
END DO
! SOLVE
CALL sls_solve(x, y, coeffs)
! FREE MEM
DEALLOCATE(x)
DEALLOCATE(x_coeffs)
DEALLOCATE(y)
!############################################################
END SUBROUTINE spline_ab_poly 



SUBROUTINE spline_ab_exp(node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: node_a
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: node_b
REAL(kind=DoubleReal), DIMENSION(:), INTENT(INOUT) :: coeffs
! Private variables
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:,:) :: x
REAL(kind=DoubleReal), ALLOCATABLE, DIMENSION(:) :: y
  
IF(SIZE(node_a,1) .EQ. 3)THEN
  ! Allocate Memory
  ALLOCATE(x(1:4,1:4))
  ALLOCATE(y(1:4))  
  x(1,1) = 1.0D0
  x(1,2) = node_a(1)
  x(1,3) = node_a(1)**2
  x(1,4) = node_a(1)**3
  x(2,1) = 1.0D0
  x(2,2) = node_b(1)
  x(2,3) = node_b(1)**2
  x(2,4) = node_b(1)**3
  x(3,1) = 0.0D0
  x(3,2) = 1.0D0
  x(3,3) = 2.0D0*node_a(1)
  x(3,4) = 3.0D0*node_a(1)**2
  x(4,1) = 0.0D0
  x(4,2) = 1.0D0
  x(4,3) = 2.0D0*node_b(1)
  x(4,4) = 3.0D0*node_b(1)**2
! make y matrix
  y(1) = log(node_a(2))
  y(2) = log(node_b(2))
  y(3) = node_a(3)/node_a(2)
  y(4) = node_b(3)/node_b(2)
ELSE IF(SIZE(node_a,1) .EQ. 4)THEN
  ! Allocate Memory
  ALLOCATE(x(1:6,1:6))
  ALLOCATE(y(1:6))  
  x(1,1) = 1.0D0
  x(1,2) = node_a(1)
  x(1,3) = node_a(1)**2
  x(1,4) = node_a(1)**3
  x(1,5) = node_a(1)**4
  x(1,6) = node_a(1)**5
  x(2,1) = 1.0D0
  x(2,2) = node_b(1)
  x(2,3) = node_b(1)**2
  x(2,4) = node_b(1)**3
  x(2,5) = node_b(1)**4
  x(2,6) = node_b(1)**5
  x(3,1) = 0.0D0
  x(3,2) = 1.0D0
  x(3,3) = 2.0D0*node_a(1)
  x(3,4) = 3.0D0*node_a(1)**2
  x(3,5) = 4.0D0*node_a(1)**3
  x(3,6) = 5.0D0*node_a(1)**4
  x(4,1) = 0.0D0
  x(4,2) = 1.0D0
  x(4,3) = 2.0D0*node_b(1)
  x(4,4) = 3.0D0*node_b(1)**2
  x(4,5) = 4.0D0*node_b(1)**3
  x(4,6) = 5.0D0*node_b(1)**4
  x(5,1) = 0.0D0
  x(5,2) = 1.0D0*node_a(3)
  x(5,3) = 2.0D0*(node_a(1)+node_a(1)*node_a(3))
  x(5,4) = 3.0D0*node_a(1)*(2.0D0*node_a(1)+node_a(1)*node_a(3))
  x(5,5) = 4.0D0*(node_a(1)**2)*(3.0D0*node_a(1)+node_a(1)*node_a(3))
  x(5,6) = 5.0D0*(node_a(1)**3)*(4.0D0*node_a(1)+node_a(1)*node_a(3))
  x(6,1) = 0.0D0
  x(6,2) = 1.0D0*node_b(3)
  x(6,3) = 2.0D0*(node_b(2)+node_b(1)*node_b(3))
  x(6,4) = 3.0D0*node_b(1)*(2.0D0*node_b(2)+node_b(1)*node_b(3))
  x(6,5) = 4.0D0*(node_b(1)**2)*(3.0D0*node_b(2)+node_b(1)*node_b(3))
  x(6,6) = 5.0D0*(node_b(1)**3)*(4.0D0*node_b(2)+node_b(1)*node_b(3))
! make y matrix
  y(1) = log(node_a(1))
  y(2) = log(node_b(2))
  y(3) = node_a(3)/node_a(1)
  y(4) = node_b(3)/node_b(2)
  y(5) = node_a(4)
  y(6) = node_b(4)
END IF
! SOLVE
CALL sls_solve(x, y, coeffs)  
! FREE MEM
DEALLOCATE(x)
DEALLOCATE(y)
END SUBROUTINE spline_ab_exp


SUBROUTINE spline_ab_exp3(node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), DIMENSION(1:3), INTENT(IN) :: node_a
REAL(kind=DoubleReal), DIMENSION(1:3), INTENT(IN) :: node_b
REAL(kind=DoubleReal), DIMENSION(1:4), INTENT(INOUT) :: coeffs
! Private variables
REAL(kind=DoubleReal), DIMENSION(1:4,1:4) :: x
REAL(kind=DoubleReal), DIMENSION(1:4) :: y
! make x 
x(1,1) = 1.0D0
x(1,2) = node_a(1)
x(1,3) = node_a(1)**2
x(1,4) = node_a(1)**3
x(2,1) = 1.0D0
x(2,2) = node_b(1)
x(2,3) = node_b(1)**2
x(2,4) = node_b(1)**3
x(3,1) = 0.0D0
x(3,2) = 1.0D0
x(3,3) = 2.0D0*node_a(1)
x(3,4) = 3.0D0*node_a(1)**2
x(4,1) = 0.0D0
x(4,2) = 1.0D0
x(4,3) = 2.0D0*node_b(1)
x(4,4) = 3.0D0*node_b(1)**2
! make y matrix
y(1) = log(node_a(2))
y(2) = log(node_b(2))
y(3) = node_a(3)/node_a(2)
y(4) = node_b(3)/node_b(2)
! SOLVE
CALL sls_solve(x, y, coeffs)  
END SUBROUTINE spline_ab_exp3





SUBROUTINE spline_ab_exp5(node_a, node_b, coeffs)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), DIMENSION(1:4), INTENT(IN) :: node_a
REAL(kind=DoubleReal), DIMENSION(1:4), INTENT(IN) :: node_b
REAL(kind=DoubleReal), DIMENSION(1:6), INTENT(INOUT) :: coeffs
! Private variables
REAL(kind=DoubleReal), DIMENSION(1:6,1:6) :: x
REAL(kind=DoubleReal), DIMENSION(1:6) :: y
! make x 
x(1,1) = 1.0D0
x(1,2) = node_a(1)
x(1,3) = node_a(1)**2
x(1,4) = node_a(1)**3
x(1,5) = node_a(1)**4
x(1,6) = node_a(1)**5
x(2,1) = 1.0D0
x(2,2) = node_b(1)
x(2,3) = node_b(1)**2
x(2,4) = node_b(1)**3
x(2,5) = node_b(1)**4
x(2,6) = node_b(1)**5
x(3,1) = 0.0D0
x(3,2) = 1.0D0
x(3,3) = 2.0D0*node_a(1)
x(3,4) = 3.0D0*node_a(1)**2
x(3,5) = 4.0D0*node_a(1)**3
x(3,6) = 5.0D0*node_a(1)**4
x(4,1) = 0.0D0
x(4,2) = 1.0D0
x(4,3) = 2.0D0*node_b(1)
x(4,4) = 3.0D0*node_b(1)**2
x(4,5) = 4.0D0*node_b(1)**3
x(4,6) = 5.0D0*node_b(1)**4
x(5,1) = 0.0D0
x(5,2) = 1.0D0*node_a(3)
x(5,3) = 2.0D0*(node_a(1)+node_a(1)*node_a(3))
x(5,4) = 3.0D0*node_a(1)*(2.0D0*node_a(1)+node_a(1)*node_a(3))
x(5,5) = 4.0D0*(node_a(1)**2)*(3.0D0*node_a(1)+node_a(1)*node_a(3))
x(5,6) = 5.0D0*(node_a(1)**3)*(4.0D0*node_a(1)+node_a(1)*node_a(3))
x(6,1) = 0.0D0
x(6,2) = 1.0D0*node_b(3)
x(6,3) = 2.0D0*(node_b(2)+node_b(1)*node_b(3))
x(6,4) = 3.0D0*node_b(1)*(2.0D0*node_b(2)+node_b(1)*node_b(3))
x(6,5) = 4.0D0*(node_b(1)**2)*(3.0D0*node_b(2)+node_b(1)*node_b(3))
x(6,6) = 5.0D0*(node_b(1)**3)*(4.0D0*node_b(2)+node_b(1)*node_b(3))
! make y matrix
y(1) = log(node_a(1))
y(2) = log(node_b(2))
y(3) = node_a(3)/node_a(1)
y(4) = node_b(3)/node_b(2)
y(5) = node_a(4)
y(6) = node_b(4)
! SOLVE
CALL sls_solve(x, y, coeffs)  
END SUBROUTINE spline_ab_exp5








