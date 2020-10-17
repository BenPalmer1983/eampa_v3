!#################################################
!# 
!#



SUBROUTINE add_es_config(rcut_in, alat_in, label_in, crystal_type_in, calc_id)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: label_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: crystal_type_in
INTEGER(KIND=StandardInteger), INTENT(OUT) :: calc_id
!###########################################################
INTEGER(KIND=StandardInteger) :: primitive_labels(1:20)
REAL(kind=DoubleReal) :: primitive_coords(1:20,1:3)
INTEGER(KIND=StandardInteger) :: atoms_per_crystal
!###########################################################

! SC = 1, BCC = 2, FCC = 3, ZB = 4
!###########################################################

! SC
IF(crystal_type_in .EQ. 1)THEN
  atoms_per_crystal = 1
  primitive_labels(1) = label_in
  primitive_coords(1,1) = 0.0D0
  primitive_coords(1,2) = 0.0D0
  primitive_coords(1,3) = 0.0D0
  !CALL add_inner(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, atoms_per_crystal, &
  !               primitive_labels(1:1), primitive_coords(1:1,:), calc_id)
  CALL add_inner(rcut_in, alat_in, label_in, crystal_type_in, atoms_per_crystal, &
                 primitive_labels(1:1), primitive_coords(1:1, :), calc_id)
! BCC
ELSE IF(crystal_type_in .EQ. 2)THEN
  atoms_per_crystal = 2
  primitive_labels(1) = label_in
  primitive_labels(2) = label_in
  primitive_coords(1,1) = 0.0D0
  primitive_coords(1,2) = 0.0D0
  primitive_coords(1,3) = 0.0D0
  primitive_coords(2,1) = 0.5D0
  primitive_coords(2,2) = 0.5D0
  primitive_coords(2,3) = 0.5D0
  CALL add_inner(rcut_in, alat_in, label_in, crystal_type_in, atoms_per_crystal, &
                 primitive_labels(1:2), primitive_coords(1:2, :), calc_id)
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
  CALL add_inner(rcut_in, alat_in, label_in, crystal_type_in, atoms_per_crystal, &
                 primitive_labels(1:4), primitive_coords(1:4, :), calc_id)
! ZB
ELSE IF(crystal_type_in .EQ. 4)THEN
  atoms_per_crystal = 8
  primitive_labels(1) = label_in
  primitive_labels(2) = label_in
  primitive_labels(3) = label_in
  primitive_labels(4) = label_in
  primitive_labels(5) = label_in
  primitive_labels(6) = label_in
  primitive_labels(7) = label_in
  primitive_labels(8) = label_in
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
  primitive_coords(5,1) = 0.25D0
  primitive_coords(5,2) = 0.25D0
  primitive_coords(5,3) = 0.25D0
  primitive_coords(6,1) = 0.75D0
  primitive_coords(6,2) = 0.75D0
  primitive_coords(6,3) = 0.25D0
  primitive_coords(7,1) = 0.25D0
  primitive_coords(7,2) = 0.75D0
  primitive_coords(7,3) = 0.75D0
  primitive_coords(8,1) = 0.75D0
  primitive_coords(8,2) = 0.25D0
  primitive_coords(8,3) = 0.75D0
  !CALL add_inner(rcut_in, alat_in, uv_in, label_in, crystal_type_in, expansion_in, atoms_per_crystal, &
  !               primitive_labels(1:8), primitive_coords(1:8,1:3), calc_id)
  CALL add_inner(rcut_in, alat_in, label_in, crystal_type_in, atoms_per_crystal, &
                 primitive_labels(1:8), primitive_coords(1:8, :), calc_id)
END IF



END SUBROUTINE add_es_config




SUBROUTINE add_inner(rcut_in, alat_in, label, crystal_type, atoms_per_crystal, &
                     primitive_labels, primitive_coords, calc_id)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rcut_in
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
INTEGER(KIND=StandardInteger), INTENT(IN) :: label
INTEGER(KIND=StandardInteger), INTENT(IN) :: crystal_type
INTEGER(KIND=StandardInteger), INTENT(IN) :: atoms_per_crystal
INTEGER(KIND=StandardInteger), INTENT(IN) :: primitive_labels(:)
REAL(kind=DoubleReal), INTENT(IN) :: primitive_coords(:,:)
INTEGER(KIND=StandardInteger), INTENT(OUT) :: calc_id
!###########################################################
INTEGER(KIND=StandardInteger) :: unmodded_labels(1:1024)
REAL(kind=DoubleReal) :: unmodded_coords(1:1024,1:3)
REAL(kind=DoubleReal) :: unmodded_uv(1:3,1:3)
INTEGER(KIND=StandardInteger) :: unmodded_atom_count
INTEGER(KIND=StandardInteger) :: slab_labels(1:1024)
REAL(kind=DoubleReal) :: slab_coords(1:1024,1:3)
REAL(kind=DoubleReal) :: slab_uv(1:3,1:3)
REAL(kind=DoubleReal) :: slab_z_mult
INTEGER(KIND=StandardInteger) :: vac_labels(1:1024)
REAL(kind=DoubleReal) :: vac_coords(1:1024,1:3)
REAL(kind=DoubleReal) :: vac_uv(1:3,1:3)
!###########################################################
REAL(kind=DoubleReal) :: strain(1:3,1:3)
REAL(kind=DoubleReal) :: uv_config(1:3,1:3)
REAL(kind=DoubleReal) :: vol_cell(1:3,1:3)
REAL(kind=DoubleReal) :: vol
REAL(kind=DoubleReal) :: sigma
REAL(kind=DoubleReal) :: padding
REAL(kind=DoubleReal) :: dz
INTEGER(KIND=StandardInteger) :: expansion
INTEGER(KIND=StandardInteger) :: n, i, j, k, m, s, cn_counter, cc_start, cc_end, ec_n, rn, cc_counter
!###########################################################

! CHECK 
n = 1
calc_id = 0
DO WHILE(calc_keys_i(n, 1) .NE. -1)
  IF(calc_keys_i(n, 1) .EQ. label .AND. &
     calc_keys_i(n, 2) .EQ. crystal_type .AND. &
     calc_keys_r(n,1) .EQ. rcut_in .AND. &
     calc_keys_r(n,2) .EQ. alat_in) THEN
    
    ! LOAD ID
    calc_id = n  
    
  END IF
  n = n + 1
END DO


IF(calc_id .EQ. 0)THEN

  !# GET BPID
  calc_configs_count = calc_configs_count + 1
  calc_id = calc_configs_count
  
  !# STORE UNIQUE DETAILS - LABEL, TYPE, RCUT, ALAT
  calc_keys_i(calc_id, 1) = label
  calc_keys_i(calc_id, 2) = crystal_type
  calc_keys_i(calc_id, 3) = cc + 1
  calc_keys_r(calc_id,1) = rcut_in
  calc_keys_r(calc_id,2) = alat_in
  
  !# 
  expansion = CEILING((2.0D0 * rcut_in) / alat_in)
  known_atoms_per_crystal(calc_id) = atoms_per_crystal
  known_expansion(calc_id) = expansion
  calc_keys_i(calc_id, 5) = expansion
  
  number_count(calc_id) = atoms_per_crystal * expansion**3
  surface_area(calc_id,1) = alat_in * alat_in * expansion  !    <-- revisit this calculation


  !# CONFIG COUNTER FOR THIS CALC_ID
  cn_counter = 0


  !########################################
  !# MAKE UNMODDED
  !########################################
  
  m = 0
  DO i=1,expansion
    DO j=1,expansion
      DO k=1,expansion
        DO n=1,SIZE(primitive_labels, 1)
          m = m + 1
          unmodded_labels(m) = primitive_labels(n)
          unmodded_coords(m,1) = (i - 1 + primitive_coords(n,1)) / expansion
          unmodded_coords(m,2) = (j - 1 + primitive_coords(n,2)) / expansion
          unmodded_coords(m,3) = (k - 1 + primitive_coords(n,3)) / expansion
        END DO
      END DO
    END DO
  END DO
  unmodded_atom_count = m
  calc_keys_i(calc_id, 6) = unmodded_atom_count
  
  ! Expand input uv
  CALL make_identity(unmodded_uv)
  unmodded_uv = expansion * unmodded_uv  
  
  ! ADD
  CALL add_config(rcut_in, alat_in, unmodded_uv, unmodded_labels(1:m), unmodded_coords(1:m,1:3))  ! Increments cc by 1
  cn_counter = cn_counter + 1
  calc_configs(calc_id, cn_counter) = cc  
  
  !MAKE NL
  CALL make_nl(cc)
  
  


  !########################################
  !# MAKE SLAB
  !########################################  
  
  ! give it r_cut on each side   
  !          rcut              rcut
  !         |    ##############    |
  
  padding = rcut_in + 2.0D0
  slab_z_mult = (expansion * alat_in + 2 * padding) / (expansion * alat_in)  
  dz = padding / (expansion * alat_in + 2 * padding)
  
  !print *, expansion, slab_z_mult
  m = 0
  DO i=1,expansion
    DO j=1,expansion
      DO k=1,expansion
        DO n=1,SIZE(primitive_labels, 1)
          m = m + 1
          slab_labels(m) = primitive_labels(n)
          slab_coords(m,1) = (i - 1 + primitive_coords(n,1)) / expansion
          slab_coords(m,2) = (j - 1 + primitive_coords(n,2)) / expansion
          slab_coords(m,3) = (k - 1 + primitive_coords(n,3)) / expansion
          
          ! Adjust Z
          slab_coords(m,3) = slab_coords(m,3) / slab_z_mult + dz
          
        END DO
      END DO
    END DO
  END DO 
  calc_keys_i(calc_id, 7) = m
  
  
  
  CALL make_identity(slab_uv)
  slab_uv(3,3) = 1.0D0 * slab_z_mult
  slab_uv = expansion * slab_uv
  CALL add_config(rcut_in, alat_in, slab_uv, slab_labels(1:m), slab_coords(1:m,1:3))  ! Increments cc by 1
  cn_counter = cn_counter + 1
  !calc_configs(calc_id, cn_counter) = cc  
  CALL make_nl(cc)
  
  
  
  !DO i=1,20
  ! Expand input uv
  !CALL make_identity(slab_uv)
  !slab_uv(3,3) = (1.0D0 + 0.01D0 * (i-10)) * slab_z_mult
  !slab_uv = expansion * slab_uv
  
  ! ADD
  !CALL add_config(rcut_in, alat_in, slab_uv, slab_labels(1:m), slab_coords(1:m,1:3))  ! Increments cc by 1
  !cn_counter = cn_counter + 1
  !calc_configs(calc_id, cn_counter) = cc  
  
  !MAKE NL
  !CALL make_nl(cc)
  !END DO
  
  
  !########################################
  !# MAKE VACANCY
  !########################################  
  
  m = 0
  DO i=1,expansion
    DO j=1,expansion
      DO k=1,expansion
        DO n=1,SIZE(primitive_labels, 1)
          IF(i .EQ. 1 .AND. j .EQ. 1 .AND. k .EQ. 1 .AND. n .EQ. 1)THEN
            ! SKIP ONE (vacancy)
          ELSE
            m = m + 1
            vac_labels(m) = primitive_labels(n)
            vac_coords(m,1) = (i - 1 + primitive_coords(n,1)) / expansion
            vac_coords(m,2) = (j - 1 + primitive_coords(n,2)) / expansion
            vac_coords(m,3) = (k - 1 + primitive_coords(n,3)) / expansion
          END IF
        END DO
      END DO
    END DO
  END DO 
  
  ! Expand input uv
  CALL make_identity(vac_uv)
  vac_uv = expansion * vac_uv * ((1.0D0 *(unmodded_atom_count - 1)) / (1.0D0 * unmodded_atom_count))**(1.0D0/3.0D0)
  
  ! ADD
  CALL add_config(rcut_in, alat_in, vac_uv, vac_labels(1:m), vac_coords(1:m,1:3))  ! Increments cc by 1
  cn_counter = cn_counter + 1
  calc_configs(calc_id, cn_counter) = cc  
  
  !MAKE NL
  CALL make_nl(cc)
  
  
  ! STORE END KEY
  calc_keys_i(calc_id, 4) = cc
  
  
  
  ! SC = 1, BCC = 2, FCC = 3, ZB = 4
  !crystal_type_in
  
END IF


END SUBROUTINE add_inner



SUBROUTINE make_identity(i_uv)
REAL(kind=DoubleReal), INTENT(INOUT) :: i_uv(1:3,1:3)
i_uv = 0.0D0
i_uv(1,1) = 1.0D0
i_uv(2,2) = 1.0D0
i_uv(3,3) = 1.0D0
END SUBROUTINE make_identity




















