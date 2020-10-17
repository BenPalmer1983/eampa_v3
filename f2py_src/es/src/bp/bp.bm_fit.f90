

SUBROUTINE fit_bm(v, e, p)
!###########################################################
Real(kind=DoubleReal), INTENT(IN) :: v(:)
Real(kind=DoubleReal), INTENT(IN) :: e(:)
Real(kind=DoubleReal), INTENT(OUT) :: p(1:4)
!###########################################################
Real(kind=DoubleReal) :: c(1:3)

Real(kind=DoubleReal) :: p_start(1:3)
Real(kind=DoubleReal) :: p_t(1:3)
Real(kind=DoubleReal) :: e_trial(1:SIZE(e,1))
Real(kind=DoubleReal) :: rss
Real(kind=DoubleReal) :: rn(1:3)
Real(kind=DoubleReal) :: J(1:SIZE(v,1),1:3)
Real(kind=DoubleReal) :: JT(1:3, 1:SIZE(v,1))
Real(kind=DoubleReal) :: JTR(1:3)
Real(kind=DoubleReal) :: H(1:3,1:3)
Real(kind=DoubleReal) :: r(1:SIZE(v,1))
Real(kind=DoubleReal) :: dp(1:3)
Real(kind=DoubleReal) :: b0p
INTEGER(KIND=StandardInteger) :: n, i
Real(kind=DoubleReal) :: rss_best
!###########################################################

! ESITMATE PARAMETERS

CALL fit(v, e, 2, c)

p(1) = (-1.0D0 * c(2)) / (2.0D0 * c(3))                      ! v0
p(2) = (c(3) * p(1) * p(1)) + (c(2) * p(1)) + c(1)   ! e0
p(3) = 2.0D0 * c(3) * p(1)                             ! b0
p(4) = 2.0D0                                           ! b0p

p_t(1:3) = p(1:3)
b0p = p(4)

DO n = 1, SIZE(e,1)
  e_trial(n) = BirchMurnMin(v(n), p_t, b0p) 
END DO
rss_best = SUM((e(:) - e_trial(:))**2)

p_start(1:3) = p(1:3)
b0p = 1.0D0
DO i=1,10
  p_t(1:3) = p_start(1:3)
  DO n=1,50
    CALL make_jacobian(p_t, b0p, v, e, J, r)
    JT = TRANSPOSE(J)
    H = MATMUL(JT, J)
    JTR = MATMUL(JT, R)
    CALL solve(H, JTR, dp)
    p_t(1:3) = p_t(1:3) + dp(1:3)
  END DO
  DO n = 1, SIZE(e,1)
    e_trial(n) = BirchMurnMin(v(n), p_t, b0p) 
  END DO
  rss = SUM((e(:) - e_trial(:))**2)
  IF(rss .LT. rss_best)THEN
    p(1:3) = p_t(1:3)
    p(4) = b0p
    rss_best = rss
  END IF
  b0p = b0p + 0.5D0
END DO

END SUBROUTINE fit_bm


SUBROUTINE make_jacobian(p, b0p, x, y, j, r)
!##################################################
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: b0p
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(OUT) :: j(1:SIZE(x,1),1:SIZE(p,1))
REAL(kind=DoubleReal), INTENT(OUT) :: r(1:SIZE(x,1))
!##################################################
INTEGER(kind=StandardInteger) :: p_size, x_size
REAL(kind=DoubleReal) :: p_temp(1:SIZE(p,1))
INTEGER(kind=StandardInteger) :: n, m
REAL(kind=DoubleReal) :: h, fa, fb
!##################################################
! Estimate dx
p_size = SIZE(p, 1)
x_size = SIZE(x, 1)
h = 0.00001D0
! Calc Jacobian
DO n=1, x_size
  fa = BirchMurnMin(x(n), p, b0p)  
  ! Store Residuals at the same time
  r(n) = y(n) - fa   
  DO m=1, p_size
    p_temp(1:p_size) = p(1:p_size)
    p_temp(m) = p_temp(m) + h    
    fb = BirchMurnMin(x(n), p_temp, b0p)
    j(n,m) = (fb-fa)/h
  END DO  
END DO

END SUBROUTINE make_jacobian   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


Function BirchMurn(volume, p) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
Implicit None  !Force declaration of all variables
! Declare variables
Real(kind=DoubleReal) :: volume
Real(kind=DoubleReal) :: p(1:4)
Real(kind=DoubleReal) :: E0, V0, B0, B0P
Real(kind=DoubleReal) :: energy
Real(kind=DoubleReal) :: eta

v0 = p(1)
e0 = p(2)
b0 = p(3)
b0p = p(4)


! Calculate energy
eta = ((1.0D0*volume)/(1.0D0*V0))**(1.0D0/3.0D0)
energy = E0+(9.0D0/16.0D0)*(B0*V0)*&
         ((eta**2-1.0D0)**2)*(6.0D0+B0P*(eta**2-1.0D0)-4.0D0*eta**2)
! Rearranged:
! energy = E+(9.0D0/16.0D0)*(B)*(&
!  V**(-1.0D0)*volume**2.0D0*(BP-4.0D0)&
!  +V**(-1.0D0/3.0D0)*volume**(4.0D0/3.0D0)*(14.0D0-3.0D0*BP)&
!  +V**(1.0D0/3.0D0)*volume**(2.0D0/3.0D0)*(3.0D0*BP-16.0D0)&
!  +V*(6.0D0-BP))
End Function BirchMurn




Function BirchMurnMin(volume, p, b0p) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
Implicit None  !Force declaration of all variables
! Declare variables
Real(kind=DoubleReal) :: volume
Real(kind=DoubleReal) :: p(1:3)
Real(kind=DoubleReal) :: E0, V0, B0, B0P
Real(kind=DoubleReal) :: energy
Real(kind=DoubleReal) :: eta

v0 = p(1)
e0 = p(2)
b0 = p(3)


! Calculate energy
eta = ((1.0D0*volume)/(1.0D0*V0))**(1.0D0/3.0D0)
energy = E0+(9.0D0/16.0D0)*(B0*V0)*&
         ((eta**2-1.0D0)**2)*(6.0D0+B0P*(eta**2-1.0D0)-4.0D0*eta**2)
! Rearranged:
! energy = E+(9.0D0/16.0D0)*(B)*(&
!  V**(-1.0D0)*volume**2.0D0*(BP-4.0D0)&
!  +V**(-1.0D0/3.0D0)*volume**(4.0D0/3.0D0)*(14.0D0-3.0D0*BP)&
!  +V**(1.0D0/3.0D0)*volume**(2.0D0/3.0D0)*(3.0D0*BP-16.0D0)&
!  +V*(6.0D0-BP))
End Function BirchMurnMin





