!#
!#  Fills array with new x,y points, also dy/dx and higher if needed
!#

SUBROUTINE speed_test()
!############################################################
REAL(kind=DoubleReal) :: d_in(1:1001,1:4)
REAL(kind=DoubleReal) :: d(1:1001,1:4)
INTEGER(kind=StandardInteger) :: n
REAL(kind=DoubleReal) :: r, yi
REAL(kind=DoubleReal) :: start_time, end_time
!############################################################


!ALLOCATE(d_in(1:1001, 1:4))
!ALLOCATE(d(1:1001, 1:4))


DO n=1,1001
  d_in(n,1) = ((n - 1.0D0) / 1000D0) * 6.0D0
  d_in(n,2) = 2.0D0 * exp(0.1D0 * d_in(n,1))
END DO

CALL fill(d_in(:,1), d_in(:,2), 1001, 4, 4, d)



CALL cpu_time(start_time)

r = 0.8D0
CALL RANDOM_NUMBER(r)
DO n=1,100000000
yi = EXP(0.2D0 * r)
CALL interpolate(r, d(:,1), d(:,2), 4, yi)
END DO

CALL cpu_time(end_time)

print *, (end_time - start_time)

END SUBROUTINE speed_test