

SUBROUTINE spline_it(spline_type, node_count, spline_x_min_in, spline_x_max_in, arr_inout)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: spline_type
INTEGER(kind=StandardInteger), INTENT(IN) :: node_count
REAL(kind=DoubleReal), INTENT(IN) :: spline_x_min_in
REAL(kind=DoubleReal), INTENT(IN) :: spline_x_max_in
REAL(kind=DoubleReal), INTENT(INOUT) :: arr_inout(:,:)
!############################################################
INTEGER(kind=StandardInteger) :: n, nn, n_start, n_end
REAL(kind=DoubleReal) :: x(1:SIZE(arr_inout,1))
REAL(kind=DoubleReal) :: y(1:SIZE(arr_inout,1))
REAL(kind=DoubleReal) :: arr_temp(1:SIZE(arr_inout,1), 1:4)
REAL(kind=DoubleReal) :: spline_x_min
REAL(kind=DoubleReal) :: spline_x_max
INTEGER(kind=StandardInteger) :: node
REAL(kind=DoubleReal) :: nx
REAL(kind=DoubleReal) :: nodes(1:node_count, 1:4)
REAL(kind=DoubleReal) :: row(1:4)
REAL(kind=DoubleReal) :: coeffs(1:node_count-1, 1:10)
INTEGER(kind=StandardInteger) :: w, l

!############################################################

l = SIZE(arr_inout,1)
w = SIZE(arr_inout,2)

x(:) = arr_inout(:,1)
y(:) = arr_inout(:,2)


spline_x_min = spline_x_min_in
IF(spline_x_min .LT. x(1))THEN
  spline_x_min = x(1)
END IF
spline_x_max = spline_x_max_in
IF(spline_x_max .GT. x(SIZE(x,1)))THEN
  spline_x_max = x(SIZE(x,1))
END IF

CALL fill(x, y, l, 4, 4, arr_temp)


DO node = 1, node_count
  nodes(node, 1) = spline_x_min + (node - 1) * ((spline_x_max - spline_x_min) / (node_count - 1))  
END DO



n_start = 0
n_end = 0
node = 1
DO n=1, SIZE(x)-1
  IF(n_start .EQ. 0 .AND. x(n) .GE. nodes(1, 1))THEN
    n_start = n
  END IF
  IF(x(n) .LE. nodes(node_count, 1))THEN
    n_end = n
  END IF

  IF(nodes(node, 1) .GE. x(n) .AND. nodes(node, 1) .LE. x(n+1))THEN    
    nn = n - 2
    IF(nn .LT. 1)THEN
      nn = 1
    ELSE IF((nn + 3) .GT. SIZE(x,1))THEN
      nn = SIZE(x,1) - 3
    END IF
    
    CALL interp4(nodes(node, 1), x(nn:nn+3), y(nn:nn+3), nodes(node, 2))
    CALL interp4dydxn(nodes(node, 1), x(nn:nn+3), y(nn:nn+3), 1, nodes(node, 3))
    CALL interp4dydxn(nodes(node, 1), x(nn:nn+3), y(nn:nn+3), 2, nodes(node, 4))
    
    node = node + 1
  END IF
END DO



  
! ADD NEW FUNCTIONS AFTER THIS POINT
  
  
END SUBROUTINE spline_it
