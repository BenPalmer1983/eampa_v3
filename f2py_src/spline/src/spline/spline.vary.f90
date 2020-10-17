SUBROUTINE vary(x, y, yvar, arr_out)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(IN) :: yvar(:)
REAL(kind=DoubleReal), INTENT(OUT) :: arr_out(1:SIZE(x,1),1:4)
!############################################################
INTEGER(kind=StandardInteger) :: n, m, i, j, nn, node_count
LOGICAL :: loop
REAL(kind=DoubleReal) :: pt(1:5*(SIZE(yvar,1)-1)+1,1:4)
REAL(kind=DoubleReal) :: nodes(1:SIZE(yvar,1),1:4)
REAL(kind=DoubleReal) :: dy, dx, c, xm
REAL(kind=DoubleReal) :: coeffs(SIZE(yvar,1)-1, 1:6)
REAL(kind=DoubleReal) :: row(1:4)
!############################################################

!# CREATE AN ARRAY OF NODES PLUS SOME INBETWEEN NODES
node_count = 5*(SIZE(yvar,1)-1)+1 
DO n=1, node_count
  pt(n, 1) = x(1) + (n - 1) * ((x(SIZE(x,1)) - x(1)) / (node_count - 1))  
END DO


i = 1
DO n=1, node_count  
  ! FIND i
  loop = .TRUE.
  IF((x(i) .LE. pt(n, 1) .AND. x(i+1) .GE. pt(n, 1)))THEN
    loop = .FALSE.
  END IF
  DO WHILE(loop)
    i = i + 1
    IF(n .GT. SIZE(x, 1))THEN
      i = SIZE(x, 1)
      loop = .FALSE.
    END IF
    IF((x(i) .LE. pt(n, 1) .AND. x(i+1) .GE. pt(n, 1)))THEN
      loop = .FALSE.
    END IF
  END DO  
  nn = i - 2
  IF((nn + 4 - 1) .GT. SIZE(x,1))THEN
    nn = SIZE(x,1) - 4 + 1
  ELSE IF(nn .LT. 1)THEN
    nn = 1
  END IF  
  CALL interp4(pt(n, 1), x(nn:nn+4-1), y(nn:nn+4-1),  pt(n, 2))    
END DO


!# VARY POSITIONS
i = 1
DO n = 1, SIZE(yvar, 1) - 1  
  dy = yvar(n+1) - yvar(n)
  dx = pt(i+4, 1) -  pt(i, 1) 
  xm = pt(i, 1) 
  DO m = 1, 5
    pt(i, 2) = pt(i, 2) +  yvar(n) + dy * ((pt(i, 1) -xm) / dx) 
    i = i + 1  
  END DO
END DO
pt(i, 2) = pt(i, 2) +  yvar(SIZE(yvar, 1))


!# FILL IN NODES
DO n=1, node_count
  nn = n - 2
  IF((nn + 4 - 1) .GT. SIZE(pt,1))THEN
    nn = SIZE(pt,1) - 4 + 1
  ELSE IF(nn .LT. 1)THEN
    nn = 1
  END IF   
  CALL interp4dydx(pt(n, 1), pt(nn:nn+4-1,1), pt(nn:nn+4-1,2),  pt(n, 3))
  CALL interp4dydxn(pt(n, 1), pt(nn:nn+4-1,1), pt(nn:nn+4-1,2), 2, pt(n, 4))
  !print *, pt(n, 1), pt(n, 2), pt(n, 3), pt(n, 4)
END DO


!# PICK OUT NODES
node_count = SIZE(yvar,1)
n = 1
DO WHILE(n .LE. node_count)
  nodes(n,:) = pt(5 * (n-1) + 1, :)
  !print *, n, nodes(n,1), nodes(n,2), nodes(n,3), nodes(n,4)
  n = n + 1
END DO

! COEFFS
DO n = 1, node_count-1
  CALL spline_ab_poly(nodes(n,:) , nodes(n+1,:) , coeffs(n, :))
  !print *,nodes(n,:) 
  !print *, nodes(n+1,:)
  !print *, coeffs(n, :)
  !print *,""
END DO

n = 1
DO i=1, SIZE(x, 1)  
  ! STORE X
  arr_out(i, 1) = x(i)
  
  ! FIND i
  loop = .TRUE.
  IF(arr_out(i, 1) .GE. nodes(n,1) .AND. arr_out(i, 1) .LE. nodes(n+1,1))THEN
    loop = .FALSE.
  END IF
  DO WHILE(loop)
    n = n + 1
    IF(i .GT. SIZE(x, 1))THEN
      n = SIZE(x, 1)
      loop = .FALSE.
    END IF
    IF(arr_out(i, 1) .GE. nodes(n,1) .AND. arr_out(i, 1) .LE. nodes(n+1,1))THEN
      loop = .FALSE.
    END IF
  END DO  
  CALL poly_row(coeffs(n, :), arr_out(i, 1), row)
  arr_out(i, 2) = row(2)
END DO

! FILL IN ARRAY
DO j = 3, 4
  DO i = 1, SIZE(arr_out, 1)  
    nn = i - 4 / 2
    IF((nn + 4 - 1) .GT. SIZE(x,1))THEN
      nn = SIZE(x,1) - 4 + 1
    ELSE IF(nn .LT. 1)THEN
      nn = 1
    END IF
    CALL interp4dydx(arr_out(i,1), arr_out(nn:nn+4-1,1), arr_out(nn:nn+4-1,j-1), arr_out(i,j))    
  END DO
END DO

END SUBROUTINE vary