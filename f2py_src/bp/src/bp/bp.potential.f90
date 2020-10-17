SUBROUTINE clear_potentials()
!###########################################################
pc = 0
mapkeys = 0
map_pn = 0
r_size = 1001
r_width = 4
pkey_temp = 0
pkey = 0
pot_rcut = 0
pot_temp = 0
pot = 0
pot_min_max = -1.0D0
potmap = 0              ! TYPE,  A,  B/Group
fgroups_dens = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
fgroups_embe = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
group_max = 0
zbl_counter = 0
zbl_i = 0            ! 1 id_1  2 id_2  3 spline_type
zbl_r = 0            ! 1 z1  2 z2  3 ra  4 rb
zbl_l = .TRUE.       ! on/off
END SUBROUTINE clear_potentials



SUBROUTINE add_potential(ftype, option_a, option_b, rcut_in, tab)
!###########################################################
INTEGER(kind=StandardInteger) :: ftype
INTEGER(kind=StandardInteger) :: option_a
INTEGER(kind=StandardInteger) :: option_b
REAL(kind=DoubleReal) :: rcut_in
REAL(kind=DoubleReal) :: tab(:,:)
!###########################################################
!###########################################################
!print *, ftype, option_a, option_b

IF(pc .EQ. 0)THEN
  r_size = SIZE(tab, 1)
END IF

pc = pc + 1

IF(pc .EQ. 1)THEN
  pkey_temp(pc, 1) = 1
ELSE
  pkey_temp(pc, 1) = pkey_temp(pc-1, 2) + 1
END IF
pkey_temp(pc, 2) = pkey_temp(pc, 1) + SIZE(tab, 1) - 1
pkey_temp(pc, 3) = ftype        ! 1 PAIR 2 DENS 3 EMBE
pkey_temp(pc, 4) = option_a     ! Type a
pkey_temp(pc, 5) = option_b     ! Type b or f group

! STORE RCUT
pot_rcut(pc) = rcut_in

! STORE POT
pot_temp(pkey_temp(pc, 1):pkey_temp(pc, 2), 1:4) = tab(:,1:4)


END SUBROUTINE add_potential




SUBROUTINE add_zbl(id1, id2, on, z1, z2, ra, rb, spline_type)
!###########################################################
INTEGER(kind=StandardInteger) :: id1
INTEGER(kind=StandardInteger) :: id2
LOGICAL :: on
REAL(kind=DoubleReal) :: z1
REAL(kind=DoubleReal) :: z2
REAL(kind=DoubleReal) :: ra
REAL(kind=DoubleReal) :: rb
INTEGER(kind=StandardInteger) :: spline_type
!###########################################################
!###########################################################

zbl_counter = zbl_counter + 1

zbl_l(zbl_counter) = on
zbl_i(zbl_counter, 1) = id1
zbl_i(zbl_counter, 2) = id2
zbl_i(zbl_counter, 3) = spline_type
zbl_r(zbl_counter, 1) = z1
zbl_r(zbl_counter, 2) = z2
zbl_r(zbl_counter, 3) = ra
zbl_r(zbl_counter, 4) = rb


END SUBROUTINE add_zbl





SUBROUTINE set_potentials()
!###########################################################
INTEGER(kind=StandardInteger) :: pn, fn, n, zn
INTEGER(kind=StandardInteger) :: mapkey, m
INTEGER(kind=StandardInteger) :: map_pn
INTEGER(kind=StandardInteger) :: a, b
REAL(kind=DoubleReal) :: rmin, rmax
LOGICAL :: loop
REAL(kind=DoubleReal) :: zbl_ra(1:4)
REAL(kind=DoubleReal) :: zbl_rb(1:4)
REAL(kind=DoubleReal) :: za, zb, ra, rb
REAL(kind=DoubleReal) :: test_fx
!###########################################################
INTEGER(kind=StandardInteger) :: fgroups(1:1000) = 0
INTEGER(kind=StandardInteger) :: x_set(1:1000) = 0
REAL(kind=DoubleReal) :: pot_arr(1:10000,1:10) = 0.0D0
INTEGER(kind=StandardInteger) :: k = 0
REAL(kind=DoubleReal) :: pot_t(1:1001,1:4) = 0.0D0
REAL(kind=DoubleReal) :: coeffs(1:10) = 0.0D0
!###########################################################


zbl_ra(1:4) = 0.0D0
zbl_rb(1:4) = 0.0D0
rmin = 0.0D0
rmax = 0.0D0
fgroups(1:1000) = 0
x_set(1:1000) = 0
pot_arr(1:10000,1:10) = 0.0D0
k = 0
pot_t(1:1001,1:4) = 0.0D0
coeffs(1:10) = 0.0D0

! REARRANGE FGROUP (IF NEEDED)
group_max = 0
DO pn =1,pc  
  IF(pkey_temp(pn, 3) .GT. 1)THEN
    fn = 0
    loop = .TRUE.
    DO WHILE(loop)
      fn = fn + 1
      IF(fgroups(fn) .EQ. pkey_temp(pn, 5) .OR. fgroups(fn) .EQ. 0)THEN
        loop = .FALSE.
      END IF
    END DO
    fgroups(fn) = pkey_temp(pn, 5)
    pkey_temp(pn, 5) = fn
    group_max = fn
  END IF
  !print *, pkey_temp(pn, 1), pkey_temp(pn, 2), pkey_temp(pn, 3), pkey_temp(pn, 4), pkey_temp(pn, 5)
END DO
!print *, group_max


! READ POTENTIALS INTO A TEMPORARY ARRAY pkey_temp
! STORE THE MINIMUM AND MAXIMUM r FOR EACH FUNCTION
! MAKE potmap(type, a, b or f) = key

pot = 0.0D0
map_pn = 0
DO pn =1,pc
  rmin = minval(pot_temp(pkey_temp(pn, 1):pkey_temp(pn, 2), 1))
  rmax = maxval(pot_temp(pkey_temp(pn, 1):pkey_temp(pn, 2), 1))  
  IF(potmap(pkey_temp(pn, 3), pkey_temp(pn, 4), pkey_temp(pn, 5)) .EQ. 0)THEN
    map_pn = map_pn + 1
    potmap(pkey_temp(pn, 3), pkey_temp(pn, 4), pkey_temp(pn, 5)) = map_pn    
    pkey(map_pn, 1) = 1 + (map_pn - 1) * r_size
    pkey(map_pn, 2) = map_pn * r_size
    pkey(map_pn, 3) = pkey_temp(pn, 3)  ! Type
    pkey(map_pn, 4) = pkey_temp(pn, 4)  ! A
    pkey(map_pn, 5) = pkey_temp(pn, 5)  ! B or F
  END IF
  mapkey = potmap(pkey_temp(pn, 3), pkey_temp(pn, 4), pkey_temp(pn, 5))
  IF(pot_min_max(mapkey, 1) .LT. 0.0D0 .OR. rmin .LT. pot_min_max(mapkey, 1))THEN
    pot_min_max(mapkey, 1) = rmin
  END IF
  IF(pot_min_max(mapkey, 2) .LT. 0.0D0 .OR. rmax .GT. pot_min_max(mapkey, 2))THEN
    pot_min_max(mapkey, 2) = rmax
  END IF 
END DO

! ADD POTENTIAL FUNCTIONS TOGETHER WITH SAME KEY (i.e. same type, a, b or f)
! INTERPOLATE SO MERGED FUNCTIONS RUN FOR SAME LENGTH rmin to rmax

DO pn =1, pc 
  mapkey = potmap(pkey_temp(pn, 3), pkey_temp(pn, 4), pkey_temp(pn, 5))
  IF(x_set(mapkey) .EQ. 0)THEN
    a = 1 + (mapkey - 1) * r_size
    b = mapkey * r_size
  END IF
  CALL fill_zoor(&
                 pot_temp(pkey_temp(pn, 1):pkey_temp(pn, 2),1), &
                 pot_temp(pkey_temp(pn, 1):pkey_temp(pn, 2),2), &
                 r_size, r_width, pot_min_max(mapkey, 1), pot_min_max(mapkey, 2), &
                 4, pot_arr(1:r_size, 1:4))
  IF(x_set(mapkey) .EQ. 0)THEN              
    x_set(mapkey) = 1
    pot(a:b, 1) = pot_arr(1:r_size, 1)
  END IF
  pot(a:b, 2:4) = pot(a:b, 2:4) +  pot_arr(1:r_size, 2:4) 
END DO
pc = MAXVAL(potmap)




! STORE F GROUPS DENS FOR EACH TYPE
DO pn =1, pc 
  IF(pkey(pn, 3) .EQ. 2)THEN
    k = 1
    DO WHILE(fgroups_dens(pkey(pn, 4), k) .NE. 0)
      k = k + 1
    END DO
    fgroups_dens(pkey(pn, 4), k) = pkey(pn, 5)    
  ELSE IF(pkey(pn, 3) .EQ. 3)THEN
    k = 1
    DO WHILE(fgroups_embe(pkey(pn, 4), k) .NE. 0)
      k = k + 1
    END DO
    fgroups_embe(pkey(pn, 4), k) = pkey(pn, 5)  
  END IF
END DO






!########################################
! APPLY ZBL

IF(zbl_counter .GT. 0)THEN
  DO pn =1, pc 
    IF(pkey(pn, 3) .EQ. 1)THEN      
      zn = 0
      DO WHILE(zn .LT. zbl_counter)
        zn = zn + 1
        IF(((zbl_i(zn, 1) .EQ. pkey(pn, 4) .AND. zbl_i(zn, 2) .EQ. pkey(pn, 5)) .OR. &
           (zbl_i(zn, 1) .EQ. pkey(pn, 5) .AND. zbl_i(zn, 2) .EQ. pkey(pn, 4))) .AND. &
           zbl_l(zn))THEN
          
          za = zbl_r(zn, 1)
          zb = zbl_r(zn, 2)
          ra = zbl_r(zn, 3)
          rb = zbl_r(zn, 4)
          
          a = pkey(pn, 1)
          b = pkey(pn, 2)
          
          ! Record "A" node of spline
          zbl_ra(1) = ra 
          CALL zblfull (za, zb, zbl_ra(1), zbl_ra(2), zbl_ra(3), zbl_ra(4))
          
          
          ! Record "B" node of spline
          zbl_rb(1) = rb          
          CALL pot_search_interpolate_a(rb, pot(a:b,1), pot(a:b,2), zbl_rb(2)) 
          CALL pot_search_interpolate_a(rb, pot(a:b,1), pot(a:b,3), zbl_rb(3)) 
          CALL pot_search_interpolate_a(rb, pot(a:b,1), pot(a:b,4), zbl_rb(4)) 
          
          ! REPLACE 0 to ra with ZBL
          loop = .TRUE.
          n = a
          DO WHILE(loop)     
            IF(n .GT. b .OR. pot(n,1) .GT. ra)THEN
              loop = .FALSE.
            ELSE
              CALL zblfull (za, zb, pot(n,1), pot(n,2), pot(n,3), pot(n,4))
            END IF
            n = n + 1
          END DO
          
          !print *, zbl_i(zn, 3)
          
          ! POLY3
          IF(zbl_i(zn, 3) .EQ. 1)THEN          
            ! SPLINE ra TO rb
            CALL spline_ab_poly(zbl_ra(1:3), zbl_rb(1:3), coeffs(1:4))          
            loop = .TRUE.
            n = a
            DO WHILE(loop)     
              IF(n .GT. b .OR. pot(n,1) .GT. rb)THEN
                loop = .FALSE.
              ELSE IF(pot(n,1) .GE. ra)THEN
                CALL poly_calc(coeffs(1:4), pot(n,1), pot(n,2))             
              END IF
              n = n + 1
            END DO
            
          ! POLY 5  
          ELSE IF(zbl_i(zn, 3) .EQ. 2)THEN          
            ! SPLINE ra TO rb
            CALL spline_ab_poly(zbl_ra(1:4), zbl_rb(1:4), coeffs(1:6))          
            loop = .TRUE.
            n = a
            DO WHILE(loop)     
              IF(n .GT. b .OR. pot(n,1) .GT. rb)THEN
                loop = .FALSE.
              ELSE IF(pot(n,1) .GE. ra)THEN
                CALL poly_calc(coeffs(1:6), pot(n,1), pot(n,2))             
              END IF
              n = n + 1
            END DO
            
          ! EXP 3  
          ELSE IF(zbl_i(zn, 3) .EQ. 3)THEN          
            ! SPLINE ra TO rb
            CALL spline_ab_exp3(zbl_ra(1:3), zbl_rb(1:3), coeffs(1:4))          
            loop = .TRUE.
            n = a
            DO WHILE(loop)     
              IF(n .GT. b .OR. pot(n,1) .GT. rb)THEN
                loop = .FALSE.
              ELSE IF(pot(n,1) .GE. ra)THEN
                CALL exp_poly_calc(coeffs(1:4), pot(n,1), pot(n,2))             
              END IF
              n = n + 1
            END DO
            
          ! EXP 5  
          ELSE IF(zbl_i(zn, 3) .EQ. 4)THEN          
            ! SPLINE ra TO rb
            CALL spline_ab_exp5(zbl_ra(1:4), zbl_rb(1:4), coeffs(1:6))          
            loop = .TRUE.
            n = a
            DO WHILE(loop)     
              IF(n .GT. b .OR. pot(n,1) .GT. rb)THEN
                loop = .FALSE.
              ELSE IF(pot(n,1) .GE. ra)THEN
                CALL exp_poly_calc(coeffs(1:6), pot(n,1), pot(n,2))             
              END IF
              n = n + 1
            END DO
          
          END IF
          
          
          ! CALC f'(x), f''(x)
          loop = .TRUE.
          n = a
          DO WHILE(loop)     
            IF(n .GT. b .OR. pot(n,1) .GT. rb)THEN
              loop = .FALSE.
            ELSE IF(pot(n,1) .GE. ra)THEN
              m = n - 2
              IF(m .LT. a)THEN
                m = a
              ELSE IF((m+3) .GT. b)THEN
                m = b - 3
              END IF    
              CALL interp4dydxn(pot(n,1), pot(m:m+3,1), pot(m:m+3,2), 1, pot(n,3))
              CALL interp4dydxn(pot(n,1), pot(m:m+3,1), pot(m:m+3,2), 2, pot(n,4))              
            END IF
            n = n + 1
          END DO          
          
          ! Break out
          zn = zbl_counter
        END IF
      END DO
    END IF
  END DO
END IF



END SUBROUTINE set_potentials




SUBROUTINE pot_search_a(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: mapkey
!###########################################################

mapkey = potmap(f_type, option_a, option_b)
CALL pot_search_interpolate_a(x, &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 1), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 2), &
                       fx)  
END SUBROUTINE pot_search_a


SUBROUTINE pot_search_b(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: mapkey
!###########################################################

mapkey = potmap(f_type, option_a, option_b)
CALL pot_search_interpolate_b(x, &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 1), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 2), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 3), &
                       fx)
  
END SUBROUTINE pot_search_b




SUBROUTINE pot_search_interpolate_a(xi, x, y, output)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(OUT) :: output
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: n_interp
INTEGER(kind=StandardInteger) :: i, j, n
REAL(kind=DoubleReal) :: li
REAL(kind=DoubleReal) :: x_min, x_max, x_range
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

! INTERP POINTS
n_interp = 4

x_min = x(1)
x_max = x(SIZE(x,1))
x_range = x_max - x_min

output = 0.0D0
IF(xi .GE. x_min .AND. xi .LE. x_max)THEN

  !# Estimate Location
  n = int((xi / x_range) * SIZE(x,1))

  !# Force within range
  IF (n .LT. 1) THEN
    n = 1
  END If
  IF (n .GT. SIZE(x,1)) THEN
    n = SIZE(x,1)
  END If

  !# Find better value if possible
  loop = .TRUE.
  DO WHILE(loop)
    IF ((xi .GE. x(n)) .AND. (xi .LE. x(n+1))) THEN
      loop = .FALSE.
    ELSE IF (x(n) .GT. xi) THEN
      n = n - 1
    ELSE IF (x(n+1) .LT. xi) THEN
      n = n + 1
    END If
    IF(n .LE. 1)THEN
      n = 1
      loop = .FALSE.
    ELSE IF(n .GE. SIZE(x,1))THEN
      n = SIZE(x,1) - 1
      loop = .FALSE.
    END If
  END DO

  !# Calculate offset
  n = n - n_interp / 2
  IF((n + n_interp - 1) .GT. SIZE(x,1))THEN
    n = SIZE(x,1) - n_interp + 1
  ELSE IF(n .LT. 1)THEN
    n = 1
  END IF

  CALL interp4(xi, x(n:n+n_interp-1), y(n:n+n_interp-1), output)
END IF
END SUBROUTINE pot_search_interpolate_a


SUBROUTINE pot_search_interpolate_b(xi, x, y, dydx, output)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(IN) :: dydx(:)
REAL(kind=DoubleReal), INTENT(OUT) :: output(1:2)
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: n_interp
INTEGER(kind=StandardInteger) :: i, j, n
REAL(kind=DoubleReal) :: li
REAL(kind=DoubleReal) :: x_min, x_max, x_range
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

! INTERP POINTS
n_interp = 4

x_min = x(1)
x_max = x(SIZE(x,1))
x_range = x_max - x_min

output = 0.0D0
IF(xi .GE. x_min .AND. xi .LE. x_max)THEN

  !# Estimate Location
  n = int((xi / x_range) * SIZE(x,1))

  !# Force within range
  IF (n .LT. 1) THEN
    n = 1
  END If
  IF (n .GT. SIZE(x,1)) THEN
    n = SIZE(x,1)
  END If

  !# Find better value if possible
  loop = .TRUE.
  DO WHILE(loop)
    IF ((xi .GE. x(n)) .AND. (xi .LE. x(n+1))) THEN
      loop = .FALSE.
    ELSE IF (x(n) .GT. xi) THEN
      n = n - 1
    ELSE IF (x(n+1) .LT. xi) THEN
      n = n + 1
    END If
    IF(n .LE. 1)THEN
      n = 1
      loop = .FALSE.
    ELSE IF(n .GE. SIZE(x,1))THEN
      n = SIZE(x,1) - 1
      loop = .FALSE.
    END If
  END DO

  !# Calculate offset
  n = n - n_interp / 2
  IF((n + n_interp - 1) .GT. SIZE(x,1))THEN
    n = SIZE(x,1) - n_interp + 1
  ELSE IF(n .LT. 1)THEN
    n = 1
  END IF

  CALL interp4(xi, x(n:n+n_interp-1), y(n:n+n_interp-1), output(1))
  CALL interp4(xi, x(n:n+n_interp-1), dydx(n:n+n_interp-1), output(2))
END IF
END SUBROUTINE pot_search_interpolate_b 


SUBROUTINE interp4(xi, x, y, yi)
! Identity for square matrix
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: y(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: yi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j
REAL(kind=DoubleReal) :: li
!############################################################
yi = 0.0D0
IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
  DO i = 1, 4
    li = 1.0D0
    DO j = 1, 4
      IF(i .NE. j) THEN
        li = li * (xi - x(j)) / (x(i) - x(j))
      END IF
    END DO
    yi = yi + li * y(i)
  END DO
END IF
!############################################################
END SUBROUTINE interp4









SUBROUTINE linspace(a, b, n, arr)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a
INTEGER(kind=StandardInteger), INTENT(IN) :: b
INTEGER(kind=StandardInteger), INTENT(IN) :: n
REAL(kind=DoubleReal), INTENT(OUT) :: arr(1:n)
!###########################################################
INTEGER(kind=StandardInteger) :: i
!###########################################################

DO i = 1, n
  arr(i) = a + (i-1) * ((b - a) / (n-1))
END DO

END SUBROUTINE linspace




SUBROUTINE zblfull (qA, qB, x, y, dy, ddy)
! parametersIn(1) = x
! parametersIn(2) = qA
! parametersIn(3) = qB
! ZBL potential, separation x, charges qA and qB
IMPLICIT NONE
!#################################
REAL(kind=DoubleReal), INTENT(IN) :: x, qA, qB
REAL(kind=DoubleReal), INTENT(OUT) :: y, dy, ddy
! Vars:  Private
Real(kind=DoubleReal) :: xVal, xs
Real(kind=DoubleReal) :: termFa, termFb, termFc, termGa, termGb, termGc
! Force none infinite result for 0
If(x.eq.0.0D0)Then
  xVal = 0.00001D0
Else
  xVal = x
End If
xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
! Calculate y
termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
0.5099D0*exp((-0.9423D0/xs)*xVal)+&
0.2802D0*exp((-0.4029D0/xs)*xVal)+&
0.02817*exp((-0.2016D0/xs)*xVal)
y = termFa * termGa
! Calculate dy
termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                 !f'(x)
termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
0.5099D0*exp((-0.9423D0/xs)*xVal)+&
0.2802D0*exp((-0.4029D0/xs)*xVal)+&
0.02817*exp((-0.2016D0/xs)*xVal)
termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&            !g'(x)
(-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
(-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
(-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
dy = termFa*termGb+termFb*termGa
! Calculate ddy
termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                        !f'(x)
termFc = (1.0D0*qA*qB)*(xVal)**(-3.0D0)*(-1.0D0)*(-2.0D0)               !f''(x)
termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                             !g(x)
0.5099D0*exp((-0.9423D0/xs)*xVal)+&
0.2802D0*exp((-0.4029D0/xs)*xVal)+&
0.02817*exp((-0.2016D0/xs)*xVal)
termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g'(x)
(-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
(-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
(-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
termGc = (-3.2D0/xs)**2*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g''(x)
(-0.9423D0/xs)**2*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
(-0.4029D0/xs)**2*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
(-0.2016D0/xs)**2*0.02817*exp((-0.2016D0/xs)*xVal)
ddy = termFa*termGc+2*termFb*termGb+termFc*termGa    
END SUBROUTINE zblfull



