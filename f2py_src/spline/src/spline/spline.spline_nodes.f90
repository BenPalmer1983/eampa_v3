SUBROUTINE spline_nodes(spline_type, nodes_in, point_count, node_expand_in, points_out)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: spline_type
REAL(kind=DoubleReal), INTENT(IN) :: nodes_in(:,:)
INTEGER(kind=StandardInteger), INTENT(IN), OPTIONAL :: node_expand_in
INTEGER(kind=StandardInteger) :: point_count
REAL(kind=DoubleReal), INTENT(OUT) :: points_out(1:point_count, 4)
!REAL(kind=DoubleReal), INTENT(IN) :: node_b(:)
!REAL(kind=DoubleReal), INTENT(OUT) :: coeffs(1:SIZE(node_a) + SIZE(node_b) - 2)
!############################################################
REAL(kind=DoubleReal), ALLOCATABLE :: nodes_interp(:,:)
REAL(kind=DoubleReal) :: x_start, x_end, x, x_inc, y
INTEGER(kind=StandardInteger) :: node_count, node_expand
INTEGER(kind=StandardInteger) :: n, m, a, b, m_last
LOGICAL :: loop
REAL(kind=DoubleReal) :: coeffs(1:10)
!############################################################


! Optional
node_expand = 1
IF(PRESENT(node_expand_in))THEN
  node_expand = node_expand_in
END IF

node_count = SIZE(nodes_in, 1)

ALLOCATE(nodes_interp(1:node_count, 1:4))

nodes_interp(:,1) = nodes_in(:,1)
nodes_interp(:,2) = nodes_in(:,2)

DO n = 1, node_count
  a = n
  b = n+3
  IF(b .GT. node_count)THEN
    a = n - 3
    b = n
  END IF
  CALL interp4dydx(nodes_in(n,1), nodes_in(a:b,1), nodes_in(a:b,2), nodes_interp(n,3))
  CALL interp4dydxn(nodes_in(n,1), nodes_in(a:b,1), nodes_in(a:b,2), 2, nodes_interp(n,4))
END DO



x_start = nodes_interp(1,1)
x_end = nodes_interp(node_count,1)
  
  
m_last = 0
m = 1
DO n = 1, point_count
  x = x_start + (n - 1) * ((x_end - x_start) / (point_count-1))
  loop = .TRUE.
  DO WHILE(loop)
    IF(x .GE. nodes_interp(m, 1) .AND. x .LE. nodes_interp(m + 1, 1))THEN
      loop = .FALSE.
    ELSE IF(m .EQ. SIZE(nodes_interp,1))THEN
      loop = .FALSE.
    ELSE  
      m = m + 1
    END IF
  END DO
  IF(spline_type .EQ. 1)THEN
    IF(m .NE. m_last)THEN    
      !CALL spline_ab(spline_type, nodes_interp(m, 1:4), nodes_interp(m+1, 1:4), coeffs(1:4))
      CALL spline_ab_poly(nodes_interp(m, 1:3), nodes_interp(m+1, 1:3), coeffs(1:4))
    END IF
    points_out(n,1) = x
    points_out(n,2) = coeffs(1) + coeffs(2) * x + coeffs(3) * x**2 + coeffs(4) * x**3
    points_out(n,3) = coeffs(2) + 2.0 * coeffs(3) * x + 3 * coeffs(4) * x**2
    points_out(n,4) = 2.0 * coeffs(3) + 6 * coeffs(4) * x
  ELSE IF(spline_type .EQ. 2)THEN
    IF(m .NE. m_last)THEN    
      !CALL spline_ab(spline_type, nodes_interp(m, 1:4), nodes_interp(m+1, 1:4), coeffs(1:4))
      CALL spline_ab_poly(nodes_interp(m, 1:4), nodes_interp(m+1, 1:4), coeffs(1:4))
    END IF
    points_out(n,1) = x
    points_out(n,2) = coeffs(1) + coeffs(2) * x + coeffs(3) * x**2 + coeffs(4) * x**3 + coeffs(5) * x**4 + coeffs(6) * x**5
    points_out(n,3) = coeffs(2) + 2.0 * coeffs(3) * x + 3.0 * coeffs(4) * x**2 + 4.0 * coeffs(5) * x**3 + 5.0 * coeffs(6) * x**4
    points_out(n,4) = 2.0 * coeffs(3) + 6.0 * coeffs(4) * x + 12.0 * coeffs(5) * x**2 + 20.0 * coeffs(6) * x**3
  END IF
  m_last = m
END DO



!CALL interpolate(xi, nodes_in(:,1), nodes_in(:,2), 4, yi


!arr_w
!fill(x, y, arr_l, arr_w)

!coeffs = 0.0D0
!IF(spline_type .EQ. 1 .AND. SIZE(node_a) .EQ. 3 .AND. SIZE(node_b) .EQ. 3)THEN
!  CALL spline_ab_poly(node_a, node_b, coeffs)
!END IF


DEALLOCATE(nodes_interp)

END SUBROUTINE spline_nodes
 
