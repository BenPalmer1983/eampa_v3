!#
!#  Fills array with new x,y points, also dy/dx and higher if needed
!#

SUBROUTINE fill(x, y, arr_l, arr_w, n_interp_in, out_arr)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: arr_l
INTEGER(kind=StandardInteger), INTENT(IN) :: arr_w
INTEGER(kind=StandardInteger), INTENT(IN), OPTIONAL :: n_interp_in
REAL(kind=DoubleReal), INTENT(OUT) :: out_arr(1:arr_l, 1:arr_w)
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: n_interp
INTEGER(kind=StandardInteger) :: i, j, n, nn
REAL(kind=DoubleReal) :: li
REAL(kind=DoubleReal) :: x_min, x_max, xi
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

! INTERP POINTS
n_interp = 4
IF(PRESENT(n_interp_in) .AND. n_interp_in .NE. 0)THEN
  n_interp = n_interp_in
END IF

x_min = minval(x)
x_max = maxval(x)

n = 1
DO i = 1, arr_l
  ! X VAL
  xi = x_min + (i - 1.0D0) * ((x_max - x_min) / (arr_l - 1.0D0))
  out_arr(i,1) = xi

  ! FIND n
  loop = .TRUE.
  IF((x(n) .LE. xi .AND. x(n+1) .GE. xi))THEN
    loop = .FALSE.
  END IF
  DO WHILE(loop)
    n = n + 1
    IF(n .GT. SIZE(x, 1))THEN
      n = SIZE(x, 1)
      loop = .FALSE.
    END IF
    IF((x(n) .LE. xi .AND. x(n+1) .GE. xi))THEN
      loop = .FALSE.
    END IF
  END DO
  
  nn = n - n_interp / 2
  IF((nn + n_interp - 1) .GT. SIZE(x,1))THEN
    nn = SIZE(x,1) - n_interp + 1
  ELSE IF(nn .LT. 1)THEN
    nn = 1
  END IF
  
  CALL interpn(xi, x(n:n+n_interp-1), y(n:n+n_interp-1), out_arr(i,2))
END DO


DO j = 3, arr_w
  DO i = 1, arr_l  
    nn = i - n_interp / 2
    IF((nn + n_interp - 1) .GT. SIZE(x,1))THEN
      nn = SIZE(x,1) - n_interp + 1
    ELSE IF(nn .LT. 1)THEN
      nn = 1
    END IF
    CALL interpndydx(out_arr(i,1), out_arr(nn:nn+n_interp-1,1), out_arr(nn:nn+n_interp-1,j-1), out_arr(i,j))    
  END DO
END DO


END SUBROUTINE fill