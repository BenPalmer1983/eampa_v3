SUBROUTINE add_alat(bp_id, alat_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: alat_in
!###########################################################
known_alat(bp_id) = alat_in
known_set(bp_id, 1) = .TRUE. 
END SUBROUTINE add_alat



SUBROUTINE add_e0(bp_id, e0_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: e0_in
!###########################################################
known_e0(bp_id) = e0_in
known_set(bp_id, 2) = .TRUE. 
END SUBROUTINE add_e0



SUBROUTINE add_b0(bp_id, b0_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: b0_in
!###########################################################
known_b0(bp_id) = b0_in
known_set(bp_id, 3) = .TRUE. 
END SUBROUTINE add_b0



SUBROUTINE add_ec(bp_id, ec_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: ec_in(1:6,1:6)
!###########################################################
known_ec(bp_id,1:6,1:6) = ec_in(1:6, 1:6)
known_set(bp_id, 4) = .TRUE. 
END SUBROUTINE add_ec



SUBROUTINE add_g(bp_id, g_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: g_in
!###########################################################
known_g(bp_id) = g_in
known_set(bp_id, 5) = .TRUE. 
END SUBROUTINE add_g



SUBROUTINE add_e(bp_id, e_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: e_in
!###########################################################
known_e(bp_id) = e_in
known_set(bp_id, 6) = .TRUE. 
END SUBROUTINE add_e



SUBROUTINE add_v(bp_id, v_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: v_in
!###########################################################
known_v(bp_id) = v_in
known_set(bp_id, 7) = .TRUE. 
END SUBROUTINE add_v










SUBROUTINE add_amu_per_crystal(bp_id, amu_per_crystal_in)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: bp_id
REAL(kind=DoubleReal), INTENT(IN) :: amu_per_crystal_in
!###########################################################
known_amu_per_crystal(bp_id) = amu_per_crystal_in
known_set(bp_id, 20) = .TRUE. 
END SUBROUTINE add_amu_per_crystal




