! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
SUBROUTINE cubic_spline_zbl_2(r, p_in, p_fixed, y)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p_in(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: Q1, Q2, rc1, rc2, rk, rs, x, e, ep, zbl_f
REAL(kind=DoubleReal) :: xa, ya, ypa, yppa, xb, yb, ypb, yppb, H
REAL(kind=DoubleReal) :: zbl_ya, zbl_yb, spline_ya, spline_yb
REAL(kind=DoubleReal) :: zbl_y, spline_y, zbl_yblend, spline_yblend
INTEGER(kind=StandardInteger) :: n_p, n_pf, join_type
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
REAL(kind=DoubleReal) :: p(1:SIZE(p_in,1))
REAL(kind=DoubleReal) :: pf(1:SIZE(p_in,1))
INTEGER(kind=StandardInteger) :: psize
REAL(kind=DoubleReal) :: dh = 1.0D-7
REAL(kind=DoubleReal) :: a(1:3)
REAL(kind=DoubleReal) :: b(1:3)
REAL(kind=DoubleReal) :: ca(1:4)
REAL(kind=DoubleReal) :: cb(1:6)
!############################################################
psize = SIZE(p,1)
y = 0.0D0
IF(SIZE(p_in,1) + 2 .EQ. SIZE(p_fixed,1))THEN  
  p(:) = p_in(1:psize)
  pf(1:psize) = p_fixed(1:psize)
  Q1 = p_fixed(psize + 1)
  Q2 = p_fixed(psize + 2)
  IF(r .EQ. 0.0D0)THEN
    y = 1.0D9
  ELSE    
    y = y + f_zbl_ackland_mendelev(r, Q1, Q2) 
    y = y + f_spline(r, p, pf, 3) 
  END IF 
END IF


END SUBROUTINE cubic_spline_zbl_2



! VECTOR SUBROUTINE
SUBROUTINE cubic_spline_zbl_2_v(r, p, p_fixed, y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
! Loop through all the values in r(:), calculate and store in y(:)
DO n = 1, SIZE(r,1)
  CALL cubic_spline_zbl_2(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_zbl_2_v



FUNCTION fade(x, m, p) RESULT (y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: x
REAL(kind=DoubleReal) :: m
REAL(kind=DoubleReal) :: p
REAL(kind=DoubleReal) :: y
!############################################################
IF(x .LE. 0.0D0)THEN
  y = 0.0D0
ELSE
  y = exp(-1.0D0 / (m * x**p))
END IF
END FUNCTION fade









