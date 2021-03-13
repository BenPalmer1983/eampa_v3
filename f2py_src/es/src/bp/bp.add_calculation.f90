

SUBROUTINE add_bp_config(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, bp_id)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
REAL(kind=DoubleReal), INTENT(IN) :: uv_in(1:3,1:3)
INTEGER(KIND=StandardInteger), INTENT(IN) :: label_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: crystal_type_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: expansion_in
INTEGER(KIND=StandardInteger), INTENT(OUT) :: bp_id
!###########################################################
INTEGER(KIND=StandardInteger) :: primitive_labels(1:4)
REAL(kind=DoubleReal) :: primitive_coords(1:4,1:3)
INTEGER(KIND=StandardInteger) :: atoms_per_crystal
!###########################################################

! SC = 1, BCC = 2, FCC = 3, ZB = 4
!###########################################################

! SC
IF(crystal_type_in .EQ. 1)THEN



! BCC
ELSE IF(crystal_type_in .EQ. 2)THEN



! FCC
ELSE IF(crystal_type_in .EQ. 3)THEN

  atoms_per_crystal = 4

  primitive_labels(1) = label_in
  primitive_labels(2) = label_in
  primitive_labels(3) = label_in
  primitive_labels(4) = label_in

  primitive_coords(1,1) = 0.0D0
  primitive_coords(1,2) = 0.0D0
  primitive_coords(1,3) = 0.0D0
  primitive_coords(2,1) = 0.5D0
  primitive_coords(2,2) = 0.5D0
  primitive_coords(2,3) = 0.0D0
  primitive_coords(3,1) = 0.5D0
  primitive_coords(3,2) = 0.0D0
  primitive_coords(3,3) = 0.5D0
  primitive_coords(4,1) = 0.0D0
  primitive_coords(4,2) = 0.5D0
  primitive_coords(4,3) = 0.5D0

  CALL add_inner(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, atoms_per_crystal, &
                 primitive_labels, primitive_coords, bp_id)

! ZB
ELSE IF(crystal_type_in .EQ. 4)THEN



END IF

END SUBROUTINE add_bp_config







SUBROUTINE add_fcc(rcut_in, alat_in, label, bp_id)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: label
INTEGER(KIND=StandardInteger), INTENT(OUT) :: bp_id
!###########################################################
INTEGER(KIND=StandardInteger) :: primitive_labels(1:4)
REAL(kind=DoubleReal) :: primitive_coords(1:4,1:3)
INTEGER(KIND=StandardInteger) :: crystal_type
INTEGER(KIND=StandardInteger) :: expansion
INTEGER(KIND=StandardInteger) :: atoms_per_crystal
REAL(kind=DoubleReal) :: uv_in(1:3,1:3)
!###########################################################

! SC = 1, BCC = 2, FCC = 3, ZB = 4
crystal_type = 3
expansion = 4
atoms_per_crystal = 4

primitive_labels(1) = label
primitive_labels(2) = label
primitive_labels(3) = label
primitive_labels(4) = label

primitive_coords(1,1) = 0.0D0
primitive_coords(1,2) = 0.0D0
primitive_coords(1,3) = 0.0D0
primitive_coords(2,1) = 0.5D0
primitive_coords(2,2) = 0.5D0
primitive_coords(2,3) = 0.0D0
primitive_coords(3,1) = 0.5D0
primitive_coords(3,2) = 0.0D0
primitive_coords(3,3) = 0.5D0
primitive_coords(4,1) = 0.0D0
primitive_coords(4,2) = 0.5D0
primitive_coords(4,3) = 0.5D0

uv_in(:,:) = 0.0D0
uv_in(1,1) = 1.0D0
uv_in(2,2) = 1.0D0
uv_in(3,3) = 1.0D0


CALL add_inner(rcut_in, alat_in, uv_in, label, crystal_type, expansion, atoms_per_crystal, &
               primitive_labels, primitive_coords, bp_id)


END SUBROUTINE add_fcc



SUBROUTINE add_inner(rcut_in, alat_in, uv_in, label, crystal_type, expansion, atoms_per_crystal, &
                     primitive_labels, primitive_coords, bp_id)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
REAL(kind=DoubleReal), INTENT(IN) :: uv_in(1:3,1:3)
INTEGER(KIND=StandardInteger), INTENT(IN) :: label
INTEGER(KIND=StandardInteger), INTENT(IN) :: crystal_type
INTEGER(KIND=StandardInteger), INTENT(IN) :: expansion
INTEGER(KIND=StandardInteger), INTENT(IN) :: atoms_per_crystal
INTEGER(KIND=StandardInteger), INTENT(IN) :: primitive_labels(:)
REAL(kind=DoubleReal), INTENT(IN) :: primitive_coords(:,:)
INTEGER(KIND=StandardInteger), INTENT(OUT) :: bp_id
!###########################################################
INTEGER(KIND=StandardInteger) :: expanded_labels(1:1024)
REAL(kind=DoubleReal) :: expanded_coords(1:1024,1:3)
REAL(kind=DoubleReal) :: strain(1:3,1:3)
REAL(kind=DoubleReal) :: uv_config(1:3,1:3)
REAL(kind=DoubleReal) :: vol_cell(1:3,1:3)
REAL(kind=DoubleReal) :: vol
REAL(kind=DoubleReal) :: sigma
INTEGER(KIND=StandardInteger) :: n, i, j, k, m, cn_counter, cc_start, cc_end, ec_n, rn, cc_counter
!###########################################################

! CHECK 
n = 1
bp_id = 0
DO WHILE(bp_keys_i(n, 1) .NE. -1)
  IF(bp_keys_i(n, 1) .EQ. label .AND. &
     bp_keys_i(n, 2) .EQ. crystal_type .AND. &
     bp_keys_r(n,1) .EQ. rcut_in .AND. &
     bp_keys_r(n,2) .EQ. alat_in) THEN
    
    ! LOAD ID
    bp_id = n  
    
  END IF
  n = n + 1
END DO


IF(bp_id .EQ. 0)THEN

  ! Calculate number of eos points
  eos_points = 2 * eos_size + 1

  !# GET BPID
  bp_configs_count = bp_configs_count + 1
  bp_id = bp_configs_count
  
  !# STORE UNIQUE DETAILS - LABEL, TYPE, RCUT, ALAT
  bp_keys_i(bp_id, 1) = label
  bp_keys_i(bp_id, 2) = crystal_type
  bp_keys_r(bp_id,1) = rcut_in
  bp_keys_r(bp_id,2) = alat_in
  
  !# 
  known_atoms_per_crystal(bp_id) = atoms_per_crystal
  known_expansion(bp_id) = expansion
  
  !# CONFIG COUNTER FOR THIS BPID
  cn_counter = 0

  m = 0
  DO i=1,expansion
    DO j=1,expansion
      DO k=1,expansion
        DO n=1,SIZE(primitive_labels, 1)
          m = m + 1
          expanded_labels(m) = primitive_labels(n)
          expanded_coords(m,1) = (i - 1 + primitive_coords(n,1)) / expansion
          expanded_coords(m,2) = (j - 1 + primitive_coords(n,2)) / expansion
          expanded_coords(m,3) = (k - 1 + primitive_coords(n,3)) / expansion
        END DO
      END DO
    END DO
  END DO
  
  
  !##################################################
  ! EQUATION OF STATE
  !##################################################
  
  cc_start = cc + 1
  cc_end = cc + eos_points
    
  cc_counter = 1
  rn = 1
  bp_cn(bp_id, rn, 1) = cc_counter
  
  DO n=1,eos_points
    cc_counter = cc_counter + 1
  
    ! Expand input uv
    uv_config = expansion * uv_in  
  
    ! Calculate strain
    sigma = (-(eos_size + 1) + n) * (eos_strain / eos_size)
    strain(:,:) = 0.0D0
    strain(1,1) = 1.0D0 + sigma
    strain(2,2) = 1.0D0 + sigma
    strain(3,3) = 1.0D0 + sigma
    
    ! Apply strain
    uv_config = matmul(strain, uv_config)

    ! Calc volume
    vol_cell(:,:) = alat_in * uv_config(:,:)
    vol = TripleProduct(vol_cell(1,:), vol_cell(2,:), vol_cell(3,:))

    CALL add_config(rcut_in, alat_in, uv_config, expanded_labels(1:m), expanded_coords(1:m,1:3))  ! Increments cc by 1
    cn_counter = cn_counter + 1
    bp_configs(bp_id, cn_counter) = cc 
    CALL make_nl(cc)

    ! RECORD DATA
    calc_strains(bp_id, 1, n) = sigma
    calc_volumes(bp_id, 1, n) = vol/m 
    
    ! REALS
    bp_data_r(cc, 1) = rcut_in              ! rcut
    bp_data_r(cc, 2) = alat_in              ! alat
    bp_data_r(cc, 3) = sigma                ! strain
    bp_data_r(cc, 4) = vol                  ! vol
    bp_data_r(cc, 5) = vol/m                ! vol/atom
    !bp_data_r(cc, 11) =                    ! epa
    !bp_data_r(cc, 12) =                    ! epa_eos
    
    ! INTS
    bp_data_i(cc, 1) = n                    ! eos number
    bp_data_i(cc, 2) = cc_start             ! bp start
    bp_data_i(cc, 3) = cc_end               ! bp end
    bp_data_i(cc, 4) = m                    ! atom count
    !print *, key(cc, 1), key(cc, 2), key(cc, 3), key(cc, 4), key(cc, 5), key(cc, 6)
  END DO

  bp_cn(bp_id, rn, 2) = cc_counter - 1
  
    
  !##################################################
  ! ELASTIC CONSTANTS
  !##################################################
  

  ! Calculate number of eos points
  ec_points = ec_size
  
  DO ec_n = 1, 9
    cc_start = cc + 1
    cc_end = cc + ec_points  
    
    rn = rn + 1
    bp_cn(bp_id, rn, 1) = cc_counter  
    
    DO n=1,ec_points
    cc_counter = cc_counter + 1
         
      ! Expand input uv
      uv_config = expansion * uv_in  
  
      ! Calculate strain
      sigma = (n - 1) * (eos_strain / ec_points) 
      strain(:,:) = 0.0D0
      
      IF(ec_n .EQ. 1)THEN
        strain(1,1) = 1.0D0 + sigma
        strain(2,2) = 1.0D0
        strain(3,3) = 1.0D0
      ELSE IF (ec_n .EQ. 2)THEN
        strain(1,1) = 1.0D0
        strain(2,2) = 1.0D0 + sigma
        strain(3,3) = 1.0D0
      ELSE IF (ec_n .EQ. 3)THEN
        strain(1,1) = 1.0D0
        strain(2,2) = 1.0D0
        strain(3,3) = 1.0D0 + sigma
      ELSE IF (ec_n .EQ. 4)THEN
        strain(1,1) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,3) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,2) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      ELSE IF (ec_n .EQ. 5)THEN
        strain(1,1) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(1,3) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,1) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      ELSE IF (ec_n .EQ. 6)THEN
        strain(1,1) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = 1.0D0 / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(1,2) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,1) = sigma / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      ELSE IF (ec_n .EQ. 7)THEN
        strain(1,1) = (1.0D0 + sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = (1.0D0 - sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = (1.0D0) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      ELSE IF (ec_n .EQ. 8)THEN
        strain(1,1) = (1.0D0 + sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = (1.0D0) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = (1.0D0 - sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      ELSE IF (ec_n .EQ. 9)THEN
        strain(1,1) = (1.0D0) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(2,2) = (1.0D0 + sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
        strain(3,3) = (1.0D0 - sigma) / ((1.0D0 - sigma**2)**(1.0D0/3.0D0))
      END IF
      
    
      ! Apply strain
      uv_config = matmul(strain, uv_config)

      ! Calc volume
      vol_cell(:,:) = alat_in * uv_config(:,:)
      vol = TripleProduct(vol_cell(1,:), vol_cell(2,:), vol_cell(3,:))

      CALL add_config(rcut_in, alat_in, uv_config, expanded_labels(1:m), expanded_coords(1:m,1:3))  ! Increments cc by 1
      cn_counter = cn_counter + 1
      bp_configs(bp_id, cn_counter) = cc 
      CALL make_nl(cc)
    
      ! RECORD DATA
      calc_strains(bp_id, ec_n + 1, n) = sigma
      calc_volumes(bp_id, ec_n + 1, n) = vol/m 
    
      ! RECORD DATA
      bp_data_r(cc, 1) = rcut_in              ! rcut
      bp_data_r(cc, 2) = alat_in              ! alat
      bp_data_r(cc, 3) = sigma                ! strain
      bp_data_r(cc, 4) = vol                  ! vol
      bp_data_r(cc, 5) = vol/m                ! vol/atom
      !bp_data_r(cc, 11) =                    ! epa
      !bp_data_r(cc, 12) =                    ! epa_eos
    
      bp_data_i(cc, 1) = n                    ! eos number
      bp_data_i(cc, 2) = cc_start             ! bp start
      bp_data_i(cc, 3) = cc_end               ! bp end
      bp_data_i(cc, 4) = m                    ! atom count
      !print *, key(cc, 1), key(cc, 2), key(cc, 3), key(cc, 4), key(cc, 5), key(cc, 6)
      
    END DO
    bp_cn(bp_id, rn, 2) = cc_counter - 1
    
    
  END DO  
  
  
END IF


END SUBROUTINE add_inner
