! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r)

! SCALAR SUBROUTINE
SUBROUTINE cubic_spline_zbl(r, p_in, p_fixed, y)
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
IF(SIZE(p_in,1) + 5 .EQ. SIZE(p_fixed,1))THEN  
  p(:) = p_in(1:psize)
  pf(1:psize) = p_fixed(1:psize)
  Q1 = p_fixed(psize + 1)
  Q2 = p_fixed(psize + 2)
  rc1 = p_fixed(psize + 3)
  rc2 = p_fixed(psize + 4)
  join_type = int(p_fixed(psize + 5))
  IF(r .EQ. 0.0D0)THEN
    y = 1.0D9
  ELSE IF(r .LT. rc1)THEN
    y = f_zbl_ackland_mendelev(r, Q1, Q2) 
  ELSE IF(r .GE. rc1 .AND. r .LE. rc2)THEN
    xa = rc1
    ya = f_zbl_ackland_mendelev(xa, Q1, Q2) 
    ypa = (f_zbl_ackland_mendelev(xa+dh, Q1, Q2) - f_zbl_ackland_mendelev(xa-dh, Q1, Q2)) / (2.0D0 * dh)
    xb = rc2
    yb = f_spline(xb, p, pf, 3)
    ypb = (f_spline(xb+dh, p, pf, 3) - f_spline(xb-dh, p, pf, 3)) / (2.0D0 * dh)
    !print *, r, f_zbl_ackland_mendelev(r, Q1, Q2), f_spline(r, p, pf, 3)
    IF(join_type .EQ. 1)THEN
      ! CUBIC
      y = f_cubic_knot_ab(r, xa, ya, ypa, xb, yb, ypb)
    ELSE IF(join_type .EQ. 2)THEN
      ! QUINTIC
      yppa = (f_zbl_ackland_mendelev(xa+dh, Q1, Q2) - 2 * f_zbl_ackland_mendelev(xa, Q1, Q2) + &
              f_zbl_ackland_mendelev(xa-dh, Q1, Q2)) / (dh**2)
      yppb = (f_zbl_ackland_mendelev(xb+dh, Q1, Q2) - 2 * f_zbl_ackland_mendelev(xb, Q1, Q2) + &
              f_zbl_ackland_mendelev(xb-dh, Q1, Q2)) / (dh**2)
      y = f_quintic_knot_ab(r, xa, ya, ypa, yppa, xb, yb, ypb, yppb)
    ELSE IF(join_type .EQ. 3)THEN
      zbl_ya = f_zbl_ackland_mendelev(rc1, Q1, Q2)
      zbl_yb = f_zbl_ackland_mendelev(rc2, Q1, Q2)
      spline_ya = f_spline(rc1, p, pf, 3)
      spline_yb = f_spline(rc2, p, pf, 3)

      zbl_y = f_zbl_ackland_mendelev(r, Q1, Q2)
      spline_y = f_spline(r, p, pf, 3)

      zbl_yblend = spline_yb + zbl_y * ((zbl_ya-spline_yb)/(zbl_ya-zbl_yb))
      spline_yblend = spline_yb + (spline_y-spline_yb)*((zbl_ya-spline_yb)/(spline_ya-spline_yb))

      y = ((rc2 - r)/(rc2-rc1))* zbl_yblend + ((r - rc1)/(rc2-rc1)) * spline_yblend
      !print *, r, zbl_y, spline_y, zbl_yblend, spline_yblend, y


    ELSE IF(join_type .EQ. 90)THEN
      a(1) = xa
      a(2) = ya
      a(3) = ypa
      b(1) = xb
      b(2) = yb
      b(3) = ypb
      CALL exp_spline_fit_a(a, b, ca)
      y = exp(ca(1) + ca(2)*x + ca(3)*x**2 + ca(4)*x**3)
    ELSE IF(join_type .EQ. 91)THEN
      a(1) = xa
      a(2) = ya
      a(3) = ypa
      b(1) = xb
      b(2) = yb
      b(3) = ypb
      CALL exp_spline_fit_b(a, b, cb)
      y = cb(1) + cb(2) * exp(cb(3) + cb(4)*x + cb(5)*x**2 + cb(6)*x**3)
    END IF
  ELSE IF(r .GT. rc2)THEN
    y = f_spline(r, p, pf, 3) 
  END IF
END IF


END SUBROUTINE cubic_spline_zbl



! VECTOR SUBROUTINE
SUBROUTINE cubic_spline_zbl_v(r, p, p_fixed, y)
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
  CALL cubic_spline_zbl(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_spline_zbl_v














