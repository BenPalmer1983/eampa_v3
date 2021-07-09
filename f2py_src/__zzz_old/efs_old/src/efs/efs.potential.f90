SUBROUTINE clear_potentials()
!###########################################################
pc = 0
mapkeys = 0
map_pn = 0
r_size = 1001
r_width = 4
pkey_input = 0
pot_input = 0
pkey_combined = 0
pot_combined = 0
pkey = 0
pot = 0
pot_rcut = 0
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



SUBROUTINE add_potential(ftype, option_a, option_b, rcut_in, tab, pykey, pc_out)
!SUBROUTINE add_potential(ftype, option_a, option_b, tab, pykey, pc_out)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: ftype
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a
INTEGER(kind=StandardInteger), INTENT(IN) :: option_b
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: tab(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: pykey
INTEGER(kind=StandardInteger), INTENT(OUT) :: pc_out
!###########################################################
!###########################################################
!print *, ftype, option_a, option_b

IF(pc .EQ. 0)THEN
  r_size = SIZE(tab, 1)
END IF

pc = pc + 1
pc_out = pc
pkey_python(pc) = pykey

IF(pc .EQ. 1)THEN
  pkey_input(pc, 1) = 1
ELSE
  pkey_input(pc, 1) = pkey_input(pc-1, 2) + 1
END IF
pkey_input(pc, 2) = pkey_input(pc, 1) + SIZE(tab, 1) - 1
pkey_input(pc, 3) = ftype        ! 1 PAIR 2 DENS 3 EMBE

! IF PAIR, A = MIN(A,B)  B = MAX(A,B)
IF(ftype .EQ. 1)THEN
  IF(option_a .GT. option_b)THEN
    pkey_input(pc, 4) = option_b     ! Type a 
    pkey_input(pc, 5) = option_a     ! Type a 
  ELSE
    pkey_input(pc, 4) = option_a     ! Type a 
    pkey_input(pc, 5) = option_b     ! Type b
  END IF
ELSE
  pkey_input(pc, 4) = option_a     ! Type a 
  pkey_input(pc, 5) = option_b     ! f group
END IF

! STORE RCUT
pot_rcut(pc) = rcut_in

! STORE POT
pot_input(pkey_input(pc, 1):pkey_input(pc, 2), 1:4) = tab(:,1:4)


END SUBROUTINE add_potential



SUBROUTINE set_potentials()
!###########################################################
INTEGER(kind=StandardInteger) :: pn, pkey, dkey, ekey, pa, pb, fn, fd_end, fe_end, n
INTEGER(kind=StandardInteger) :: counter(1:max_potfunctions) = 0
INTEGER(kind=StandardInteger) :: key_list(1:max_potfunctions, 1:10) = 0
INTEGER(kind=StandardInteger) :: at, bt, ac, bc, ac_end
REAL(kind=DoubleReal) :: f_min(1:max_potfunctions) = 0.0D0
REAL(kind=DoubleReal) :: f_max(1:max_potfunctions) = 0.0D0
!###########################################################
REAL(kind=DoubleReal) :: pot_t(1:10000,1:4) = 0.0D0
!###########################################################


! ZERO ARRAYS
pot_key(:,:) = 0
pot_data(:,:) = 0.0D0


!##########
! LOOP 1
!##########

! SET label_max AND fgroup_max
! Calculate offsets dens_key_offset and embe_key_offset
fd_end = SIZE(fgroups_dens, 2)
fe_end = SIZE(fgroups_embe, 2)
DO pn =1,pc  
  ! Get label_max and fgroup_max
  IF(pkey_input(pn, 3) .EQ. 1)THEN
    label_max = MAX(label_max, pkey_input(pn, 4)) ! LABEL A
    label_max = MAX(label_max, pkey_input(pn, 5)) ! LABEL B
  ELSE
    label_max = MAX(label_max, pkey_input(pn, 4))   ! LABEL A
    fgroup_max = MAX(fgroup_max, pkey_input(pn, 5)) ! FGROUP
  END IF

  ! Store fgroups for each label type (Density)
  IF(pkey_input(pn, 3) .EQ. 2)THEN
    DO n = 1, fd_end
      IF(fgroups_dens(pkey_input(pn, 4), n) .EQ. pkey_input(pn, 5))THEN 
        EXIT 
      ELSE IF(fgroups_dens(pkey_input(pn, 4), n) .EQ. 0)THEN
        fgroups_dens(pkey_input(pn, 4), n) = pkey_input(pn, 5)
        EXIT
      END IF
    END DO
  END IF

  ! Store fgroups for each label type (Embedding)
  IF(pkey_input(pn, 3) .EQ. 3)THEN
    DO n = 1, fe_end
      IF(fgroups_embe(pkey_input(pn, 4), n) .EQ. pkey_input(pn, 5))THEN 
        EXIT 
      ELSE IF(fgroups_embe(pkey_input(pn, 4), n) .EQ. 0)THEN
        fgroups_embe(pkey_input(pn, 4), n) = pkey_input(pn, 5)
        EXIT
      END IF
    END DO
  END IF
    
END DO



CALL pair_key(label_max, label_max, dens_key_offset)
embe_key_offset = dens_key_offset + label_max * fgroup_max


! Check for potentials that need to be combined, and store min/max values
DO pn =1,pc  
  at = pkey_input(pn, 1)
  bt = pkey_input(pn, 2)
  CALL get_pot_key(pkey_input(pn, 3), pkey_input(pn, 4), pkey_input(pn, 5), pkey)
  counter(pkey) = counter(pkey) + 1
  key_list(pkey, counter(pkey)) = pn

  IF(counter(pkey) .EQ. 1)THEN
    f_min(pkey) = pot_input(at, 1)     ! function min
    f_max(pkey) = pot_input(bt, 1)     ! function max
  ELSE
    f_min(pkey) = MIN(f_min(pkey), pot_input(at, 1))
    f_max(pkey) = MAX(f_max(pkey), pot_input(bt, 1))
  END IF
END DO


ac_end = 1
DO pn =1,pc   
  CALL get_pot_key(pkey_input(pn, 3), pkey_input(pn, 4), pkey_input(pn, 5), pkey)
  ! Get start/end from input data
  at = pkey_input(pn, 1)
  bt = pkey_input(pn, 2)
  IF(pot_key(pkey, 1) .EQ. 0)THEN
    ac = ac_end
    bc = ac + r_size - 1
    ac_end = ac_end + r_size
    pot_key(pkey, 1) = ac
    pot_key(pkey, 2) = bc
  ELSE
    ac = pot_key(pkey, 1)
    bc = pot_key(pkey, 2)
  END IF
  IF(counter(pkey) .EQ. 1)THEN
    ! FILL TEMPORARY ARRAY
    pot_t(1:10000,1:4) = 0.0D0
    CALL fill(pot_input(at:bt, 1), pot_input(at:bt, 2), r_size, r_width, &
              r_interp, pot_t(1:r_size, :))      
    ! STORE IN PAIR ARRAY
    pot_data(ac:bc,1:r_width) = pot_t(1:r_size, 1:r_width)
  ELSE
    ! FILL TEMPORARY ARRAY
    pot_t(1:10000,1:4) = 0.0D0
    CALL fill_zoor(pot_input(at:bt, 1), pot_input(at:bt, 2), r_size, r_width, &
              f_min(pkey), f_max(pkey), &
              r_interp, pot_t(1:r_size, :))      
    ! ADD TO PAIR ARRAY
    pot_data(ac:bc,1:r_width) = pot_data(ac:bc,1:r_width) + pot_t(1:r_size, 1:r_width)      
  END IF  
END DO

END SUBROUTINE set_potentials






! KEYS
SUBROUTINE get_pot_key(ftype, a, b, key)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: ftype, a, b
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
!###########################################################
IF(ftype .EQ. 1)THEN
  CALL pair_key(a, b, key)
ELSE IF(ftype .EQ. 2)THEN
  CALL dens_key(a, b, key)
ELSE IF(ftype .EQ. 3)THEN
  CALL embe_key(a, b, key)
END IF

END SUBROUTINE get_pot_key





! KEYS
SUBROUTINE pair_key(a, b, key)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, b
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
! Min/Max
IF(a .GT. label_max .OR. b .GT. label_max)THEN
  key = 0
ELSE
  min_key = a
  max_key = b
  IF(a .GT. b)THEN
    max_key = a
    min_key = b
  END IF
  key = (max_key*(max_key-1))/2 + min_key
END IF
END SUBROUTINE pair_key

! KEYS
SUBROUTINE dens_key(a, fgroup, key)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, fgroup
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
IF(a .GT. label_max .OR. fgroup .GT. fgroup_max)THEN
  key = 0
ELSE
  key = dens_key_offset + fgroup + (a - 1) * fgroup_max
END IF
END SUBROUTINE dens_key

! KEYS
SUBROUTINE embe_key(a, fgroup, key)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, fgroup
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
IF(a .GT. label_max .OR. fgroup .GT. fgroup_max)THEN
  key = 0
ELSE
  key = dens_key_offset + label_max * fgroup_max + fgroup + (a - 1) * fgroup_max
END IF
END SUBROUTINE embe_key






SUBROUTINE pot_search(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(f_type, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_a(x, pot_data(a:b, 1), pot_data(a:b, 2), fx)  
END IF
END SUBROUTINE pot_search

SUBROUTINE pot_search_pair(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(1, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_a(x, pot_data(a:b, 1), pot_data(a:b, 2), fx)  
END IF
END SUBROUTINE pot_search_pair

SUBROUTINE pot_search_dens(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(2, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_a(x, pot_data(a:b, 1), pot_data(a:b, 2), fx)  
END IF
END SUBROUTINE pot_search_dens

SUBROUTINE pot_search_embe(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(3, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_a(x, pot_data(a:b, 1), pot_data(a:b, 2), fx)  
END IF
END SUBROUTINE pot_search_embe



SUBROUTINE pot_search_arr(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(f_type, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_b(x, pot_data(a:b, 1), pot_data(a:b, 2), pot_data(a:b, 3), fx)
END IF
END SUBROUTINE pot_search_arr

SUBROUTINE pot_search_pair_arr(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(1, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_b(x, pot_data(a:b, 1), pot_data(a:b, 2), pot_data(a:b, 3), fx)
END IF
END SUBROUTINE pot_search_pair_arr

SUBROUTINE pot_search_dens_arr(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(2, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_b(x, pot_data(a:b, 1), pot_data(a:b, 2), pot_data(a:b, 3), fx)
END IF
END SUBROUTINE pot_search_dens_arr

SUBROUTINE pot_search_embe_arr(option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(3, option_a, option_b, pkey)
a = pot_key(pkey, 1)
b = pot_key(pkey, 2)
IF(a .GT. 0 .AND. b .GT. a)THEN
  CALL pot_search_interpolate_b(x, pot_data(a:b, 1), pot_data(a:b, 2), pot_data(a:b, 3), fx)
END IF
END SUBROUTINE pot_search_embe_arr


















SUBROUTINE pot_search_a(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: mapkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
a = option_a
b = option_b
IF(f_type .EQ. 1)THEN
  IF(option_a .GT. option_b)THEN
    b = option_a
    a = option_b
  END IF
END IF
mapkey = potmap(f_type, a, b)
IF(mapkey .EQ. 0)THEN   ! RETURN ZERO IF NO DATA/POTENTIAL
  fx = 0.0
ELSE
  CALL pot_search_interpolate_a(x, &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 1), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 2), &
                       fx)  
END IF
END SUBROUTINE pot_search_a


SUBROUTINE pot_search_b(f_type, option_a, option_b, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx(1:2)
!###########################################################
INTEGER(kind=StandardInteger) :: mapkey
INTEGER(kind=StandardInteger) :: a, b
!###########################################################
a = option_a
b = option_b
IF(f_type .EQ. 1)THEN
  IF(option_a .GT. option_b)THEN
    b = option_a
    a = option_b
  END IF
END IF
mapkey = potmap(f_type, a, b)
IF(mapkey .EQ. 0)THEN   ! RETURN ZERO IF NO DATA/POTENTIAL
  fx = 0.0
ELSE
  CALL pot_search_interpolate_b(x, &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 1), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 2), &
                       pot(pkey(mapkey, 1):pkey(mapkey, 2), 3), &
                       fx)
END IF  
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


SUBROUTINE interp4dydx(xi, x, y, ypi)
! Interpolate and return derivative at xi
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: y(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: ypi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, k, n
REAL(kind=DoubleReal) :: fx, gx, psum
!############################################################
n = 4
IF (SIZE(x,1) .EQ. SIZE(y,1) .AND. n .EQ. SIZE(x,1)) THEN
  ypi = 0.0D0
  Do i=1,SIZE(x,1)
    fx = 1.0D0
    gx = 0.0D0
    Do j=1,SIZE(x,1)
      If(i .NE. j) Then
        fx = fx / (x(i) - x(j))
        psum = 1.0D0
        Do k=1,SIZE(x,1)
          If((i .NE. k) .AND. (j .NE. k))Then
            psum = psum * (xi - x(k))
          End If
        End Do
        gx = gx + psum
      End If
    End Do
    ypi = ypi + fx * gx * y(i)
  End Do
END IF
!############################################################
END SUBROUTINE interp4dydx






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





