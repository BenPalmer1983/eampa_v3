SUBROUTINE get_nodes(x, y, x_min, x_max, node_count, nodes)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(IN) :: x_min
REAL(kind=DoubleReal), INTENT(IN) :: x_max
INTEGER(kind=StandardInteger), INTENT(IN) :: node_count
REAL(kind=DoubleReal), INTENT(OUT) :: nodes(1:node_count,1:4)
!############################################################
INTEGER(kind=StandardInteger) :: n, i, nn
LOGICAL :: loop
!############################################################

DO n=1, node_count
  nodes(n, 1) = x_min + (n - 1) * ((x_max - x_min) / (node_count - 1))
END DO

i = 1
DO n=1, node_count  
  ! FIND i
  loop = .TRUE.
  IF((x(i) .LE. nodes(n, 1) .AND. x(i+1) .GE. nodes(n, 1)))THEN
    loop = .FALSE.
  END IF
  DO WHILE(loop)
    i = i + 1
    IF(n .GT. SIZE(x, 1))THEN
      i = SIZE(x, 1)
      loop = .FALSE.
    END IF
    IF((x(i) .LE. nodes(n, 1) .AND. x(i+1) .GE. nodes(n, 1)))THEN
      loop = .FALSE.
    END IF
  END DO
  
  nn = i - 2
  IF((nn + 4 - 1) .GT. SIZE(x,1))THEN
    nn = SIZE(x,1) - 4 + 1
  ELSE IF(nn .LT. 1)THEN
    nn = 1
  END IF
  
  CALL interp4(nodes(n, 1), x(nn:nn+4-1), y(nn:nn+4-1),  nodes(n, 2))
  CALL interp4dydx(nodes(n, 1), x(nn:nn+4-1), y(nn:nn+4-1),  nodes(n, 3))
  CALL interp4dydxn(nodes(n, 1), x(nn:nn+4-1), y(nn:nn+4-1), 2, nodes(n, 4))
    
  !print *, n, nodes(n, 1), nodes(n, 2), nodes(n, 3), nodes(n, 4)
END DO






END SUBROUTINE get_nodes


