!##########################################################
! Add config and increments config counter cc
! Makes ghost cell


SUBROUTINE add_config(rcut_in, alat_in, uv_in, labels_in, coords_in)
!###########################################################
REAL(kind=DoubleReal) :: rcut_in
REAL(kind=DoubleReal) :: alat_in
REAL(kind=DoubleReal) :: uv_in(1:3,1:3)
INTEGER(kind=StandardInteger) :: labels_in(:)
REAL(kind=DoubleReal) :: coords_in(:,:)

!###########################################################
INTEGER(kind=StandardInteger) :: ccz
INTEGER(kind=StandardInteger) :: ka, kb, kc, kd
REAL(kind=DoubleReal) :: to_cart(1:3, 1:3)

INTEGER(KIND=StandardInteger) :: cn, an, cx, cy, cz, n, m, i


REAL(kind=DoubleReal) :: shift(1:3)
REAL(kind=DoubleReal) :: coords_crystal(1:3)
REAL(kind=DoubleReal) :: coords_cartesian(1:3)
REAL(kind=DoubleReal) :: half(1:3)
REAL(kind=DoubleReal) :: c(1:3)
REAL(kind=DoubleReal) :: ct(1:3)
REAL(kind=DoubleReal) :: td(1:3)
REAL(kind=DoubleReal) :: rd(1:3) 
REAL(kind=DoubleReal) :: r(1:3) 
REAL(kind=DoubleReal) :: vol_cell(1:3,1:3)
!###########################################################

cc = cc + 1
ccz = cc - 1

rcut(cc) = rcut_in
alat(cc) = alat_in
uv(3 * ccz + 1: 3 * ccz + 3, 1:3) = uv_in(1:3,1:3) 

IF(SIZE(labels_in,1) .EQ. SIZE(coords_in,1))THEN

  ! STORE KEYS
  IF(cc .EQ. 1)THEN
    ka = 1
    kb = ka + SIZE(labels_in, 1) - 1
  ELSE  
    ka = key(cc-1, 2) + 1
    kb = ka + SIZE(labels_in, 1) - 1
  END IF
  key(cc, 1) = ka
  key(cc, 2) = kb
  
  ! STORE
  labels(ka:kb) = labels_in(:)
  coords(ka:kb, 1:3) = coords_in(:, 1:3)
      
  ! MAKE CARTESIAN 
  m = 0
  DO n = ka, kb
    m = m + 1
    ids(n) = m
    DO i = 1,3
      CALL Modulus(coords(n, i), 1.0D0, coords(n, i))
    END DO
    coords(n, 4:6) = matmul(alat(cc) * uv(3 * ccz + 1: 3 * ccz + 3, 1:3), coords(n, 1:3))
  END DO  
  
  ! MAKE GHOST 
  IF(cc .EQ. 1)THEN
    kc = 1
  ELSE  
    kc = key(cc-1, 4) + 1
  END IF  
  
  half = 0.5D0
  c = alat(cc) * MATMUL(half, uv(3 * ccz + 1: 3 * ccz + 3, 1:3))
  td = gthreshold * (c + rcut(cc))
  
  vol_cell(:,:) =  alat(cc) * uv(3 * ccz + 1: 3 * ccz + 3, 1:3)
  volume(cc) = TripleProduct(vol_cell(1,:), vol_cell(2,:), vol_cell(3,:))
  
  an = kc

  DO cx = -1, 1
    DO cy = -1 ,1
      DO cz = -1, 1
        IF(cx .EQ. 0 .AND. cy .EQ. 0 .AND. cz .EQ. 0) THEN
          DO n = ka, kb
            ghostids(an) = ids(n)
            ghostlabels(an) = labels(n)
            ghostcoords(an, 1:3) = coords(n, 4:6)
            ghosthalo(an) = .FALSE. 
            an = an + 1
          END DO
        ELSE     
          ! Shift
          shift(1) = 1.0D0 * cx
          shift(2) = 1.0D0 * cy
          shift(3) = 1.0D0 * cz
        
          DO n = ka, kb 
            ct(1:3) = matmul(alat(cc) * uv(3 * ccz + 1: 3 * ccz + 3, 1:3), coords(n, 1:3) + shift(1:3))
            rd(1:3) = abs(ct(1:3) - c(1:3))          
            IF(rd(1) .LE. td(1) .AND. rd(2) .LE. td(2) .AND. rd(3) .LE. td(3))THEN
              ghostids(an) = ids(n)
              ghostlabels(an) = labels(n)
              ghostcoords(an, 1:3) = ct(1:3)
              ghosthalo(an) = .TRUE.  
              an = an + 1
            END IF
          END DO 
        END IF   
      END DO
    END DO
  END DO
  key(cc, 3) = kc
  key(cc, 4) = an - 1  
END IF

END SUBROUTINE add_config 

