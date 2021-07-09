
!
!   Fits 
!   y = c + m exp(b0 + b1 x + b2 x**2 + b3 x**3)

SUBROUTINE exp_spline_fit_a(a_in, b_in, c_out)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a_in(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal), INTENT(IN) :: b_in(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal), INTENT(OUT) :: c_out(1:4)  ! c, m, b0, b1, b2, b3  
!############################################################
REAL(kind=DoubleReal) :: a(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal) :: b(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal) :: c(1:4)    ! b0, b1, b2, b3  
INTEGER(kind=StandardInteger) :: i, k, n
REAL(kind=DoubleReal) :: m
REAL(kind=DoubleReal) :: yc
REAL(kind=DoubleReal) :: h
REAL(kind=DoubleReal) :: J(1:4,1:4)
REAL(kind=DoubleReal) :: JT(1:4,1:4)
REAL(kind=DoubleReal) :: JTJ(1:4,1:4)
REAL(kind=DoubleReal) :: JTJ_inv(1:4,1:4)
REAL(kind=DoubleReal) :: JTR(1:4)
REAL(kind=DoubleReal) :: R(1:4)
REAL(kind=DoubleReal) :: rf(1:4)
REAL(kind=DoubleReal) :: rb(1:4)
REAL(kind=DoubleReal) :: cf(1:4)
REAL(kind=DoubleReal) :: cb(1:4)
REAL(kind=DoubleReal) :: dp(1:4)
REAL(kind=DoubleReal) :: rss_this = 0.0D0
REAL(kind=DoubleReal) :: rss_last = 0.0D0
REAL(kind=DoubleReal) :: c_last(1:4)
!############################################################

a(1:3) = a_in(1:3)
b(1:3) = b_in(1:3)

c_out(:) = 0.0D0

c(1:4) = 0.0D0

J = 0.0D0
DO n = 1,20
  R = residual(a, b, c(1:4))
  h = 1.0D-8
  DO i = 1, 4
    cf(1:4) = c(1:4)
    cf(i) = cf(i) + 0.5D0 * h
    rf = residual(a, b, cf(1:4))
	  cb(1:4) = c(1:4)
	  cb(i) = cb(i) - 0.5D0 * h
    rb = residual(a, b, cb(1:4))
    DO k =1,4
      J(k, i) = (rf(k) - rb(k)) / h
    END DO    
  END DO
	JT = TRANSPOSE(J)
	JTJ = MATMUL(JT, J)
	JTR = MATMUL(JT, R)
  CALL sls_solve(JTJ, -JTR, dp)
  c(1:4) = c(1:4) + dp(1:4)
  rss_this =  rss(a, b, c)
  IF(n .GT. 1)THEN
    IF(abs(rss_this - rss_last) < 1.0D-10)THEN 
	  	c_out(1:4) = c(1:4)
		  RETURN
    ELSE IF(rss_this > rss_last)THEN 
      c_out(1:4) = c_last(1:4)
      RETURN
    END IF
  END IF
  rss_last = rss_this
  c_last = c
END DO
c_out(1:4) = c(1:4)
RETURN
END SUBROUTINE exp_spline_fit_a


!
!   Fits 
!   y = c + m exp(b0 + b1 x + b2 x**2 + b3 x**3)

SUBROUTINE exp_spline_fit_b(a_in, b_in, c_out)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a_in(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal), INTENT(IN) :: b_in(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal), INTENT(OUT) :: c_out(1:6)  ! c, m, b0, b1, b2, b3  
!############################################################
REAL(kind=DoubleReal) :: a(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal) :: b(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal) :: c(1:4)    ! b0, b1, b2, b3  
INTEGER(kind=StandardInteger) :: i, k, n
REAL(kind=DoubleReal) :: m
REAL(kind=DoubleReal) :: yc
REAL(kind=DoubleReal) :: h
REAL(kind=DoubleReal) :: J(1:4,1:4)
REAL(kind=DoubleReal) :: JT(1:4,1:4)
REAL(kind=DoubleReal) :: JTJ(1:4,1:4)
REAL(kind=DoubleReal) :: JTJ_inv(1:4,1:4)
REAL(kind=DoubleReal) :: JTR(1:4)
REAL(kind=DoubleReal) :: R(1:4)
REAL(kind=DoubleReal) :: rf(1:4)
REAL(kind=DoubleReal) :: rb(1:4)
REAL(kind=DoubleReal) :: cf(1:4)
REAL(kind=DoubleReal) :: cb(1:4)
REAL(kind=DoubleReal) :: dp(1:4)
REAL(kind=DoubleReal) :: rss_this = 0.0D0
REAL(kind=DoubleReal) :: rss_last = 0.0D0
REAL(kind=DoubleReal) :: c_last(1:4)
!############################################################

a(1:3) = a_in(1:3)
b(1:3) = b_in(1:3)

m = a(2) - b(2)
yc = b(2)

a(2) = 1.0D0
b(2) = 0.0D0
a(3) = a(3) / abs(m)
b(3) = b(3) / abs(m)

c_out(:) = 0.0D0
c_out(1) = yc
c_out(2) = m

c(1:4) = 0.0D0

J = 0.0D0
DO n = 1,20
  R = residual(a, b, c(1:4))
  h = 1.0D-8
  DO i = 1, 4
    cf(1:4) = c(1:4)
    cf(i) = cf(i) + 0.5D0 * h
    rf = residual(a, b, cf(1:4))
	  cb(1:4) = c(1:4)
	  cb(i) = cb(i) - 0.5D0 * h
    rb = residual(a, b, cb(1:4))
    DO k =1,4
      J(k, i) = (rf(k) - rb(k)) / h
    END DO    
  END DO
	JT = TRANSPOSE(J)
	JTJ = MATMUL(JT, J)
	JTR = MATMUL(JT, R)
  CALL sls_solve(JTJ, -JTR, dp)
  c(1:4) = c(1:4) + dp(1:4)
  rss_this =  rss(a, b, c)
  IF(n .GT. 1)THEN
    IF(abs(rss_this - rss_last) < 1.0D-10)THEN 
	  	c_out(3:6) = c(1:4)
		  RETURN
    ELSE IF(rss_this > rss_last)THEN 
      c_out(3:6) = c_last(1:4)
      RETURN
    END IF
  END IF
  rss_last = rss_this
  c_last = c
END DO
c_out(3:6) = c(1:4)
RETURN
END SUBROUTINE exp_spline_fit_b





FUNCTION f_exp(x, c) RESULT (y)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: x
REAL(kind=DoubleReal) :: c(1:4)
REAL(kind=DoubleReal) :: y
!############################################################
!############################################################
y = exp(c(1) + c(2)*x + c(3)*x**2 + c(4)*x**3)
END FUNCTION f_exp


FUNCTION fp_exp(x, c) RESULT (y)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: x
REAL(kind=DoubleReal) :: c(1:4)
REAL(kind=DoubleReal) :: y
!############################################################
!############################################################
y = (c(2) + 2.0D0 * c(3)*x + 3.0D0 * c(4)*x**2) * exp(c(1) + c(2)*x + c(3)*x**2 + c(4)*x**3)
END FUNCTION fp_exp



FUNCTION f(x, c) RESULT (y)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: x
REAL(kind=DoubleReal) :: c(1:6)
REAL(kind=DoubleReal) :: y
!############################################################
!############################################################
y = c(1) + c(2) * exp(c(3) + c(4)*x + c(5)*x**2 + c(6)*x**3)
END FUNCTION f



FUNCTION fv(x, c) RESULT (y)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: x(:)
REAL(kind=DoubleReal) :: c(1:6)
REAL(kind=DoubleReal) :: y(1:SIZE(x,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n=1, SIZE(x,1)
  y(n) = f(x(n), c)
END DO
END FUNCTION fv


FUNCTION rss(a, b, c) RESULT (rss_v)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: a(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal) :: b(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal) :: c(1:4)    ! b0, b1, b2, b3  
!############################################################
REAL(kind=DoubleReal) :: rss_v
!############################################################
rss_v = 0.0D0
rss_v = rss_v + (a(2) - f_exp(a(1), c(1:4)))**2
rss_v = rss_v + (a(3) - fp_exp(a(1), c(1:4)))**2
rss_v = rss_v + (b(2) - f_exp(b(1), c(1:4)))**2
rss_v = rss_v + (b(3) - fp_exp(b(1), c(1:4)))**2
END FUNCTION rss





FUNCTION residual(a, b, c) RESULT (res)
!############################################################
! residual
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: a(1:3)    ! xa, ya, ypa  
REAL(kind=DoubleReal) :: b(1:3)    ! xb, yb, ypb  
REAL(kind=DoubleReal) :: c(1:4)    ! b0, b1, b2, b3  
!############################################################
REAL(kind=DoubleReal) :: res(1:4)
!############################################################
res(1:4) = 0.0D0
res(1) = (a(2) - f_exp(a(1), c(1:4)))**2
res(2) = (a(3) - fp_exp(a(1), c(1:4)))**2
res(3) = (b(2) - f_exp(b(1), c(1:4)))**2
res(4) = (b(3) - fp_exp(b(1), c(1:4)))**2
END FUNCTION residual
















