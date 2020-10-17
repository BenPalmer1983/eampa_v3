SUBROUTINE calculate_bp()
!###########################################################
INTEGER(KIND=StandardInteger) :: bp_id
!###########################################################

! Reset total rss over all 
rss_total_rss = 0.0D0
rss_total_rss_w = 0.0D0

! Loop through
bp_id = 1
DO WHILE(bp_keys_i(bp_id, 1) .NE. -1)
  CALL calculate_bp_inner(bp_id)
  bp_id = bp_id + 1
END DO

END SUBROUTINE calculate_bp


SUBROUTINE calculate_bp_inner(bp_id)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: bp_id
!###########################################################
INTEGER(KIND=StandardInteger) :: a, b, s, n, rn, ec_n, i, j
Real(kind=DoubleReal) :: p(1:4)
Real(kind=DoubleReal) :: c(1:3)
Real(kind=DoubleReal) :: c11, c22, c33, c44, c55, c66, c12, c13, c23
Real(kind=DoubleReal) :: h, k, pi, na, g, b0
!###########################################################

! GET KEYS
rn = 1
a = bp_cn(bp_id, rn, 1)
b = bp_cn(bp_id, rn, 2)
s = b - a + 1


! STORE ENERGIES, FIT, STORE CALC VALUES
calc_energies(bp_id, 1, 1:s) = config_energy(a:b, 6)
calc_sizes(bp_id, 1) = s
!print *, calc_energies(bp_id, 1, 1:s)
CALL fit_bm(calc_volumes(bp_id, 1, 1:s), calc_energies(bp_id, 1, 1:s), p(1:4))
DO n = 1, s
  calc_energies_fit(bp_id, 1, n) = BirchMurn(calc_volumes(bp_id, 1, n), p(1:4))  
END DO

! STORE BPS
calc_v0(bp_id) = p(1)
calc_alat(bp_id) = (known_atoms_per_crystal(bp_id) * p(1))**(1.0D0/3.0D0)
calc_e0(bp_id) = p(2)
calc_b0(bp_id) = p(3)


! DENSITY
calc_rho(bp_id) = (1660.54D0 * known_amu_per_crystal(bp_id)) / (calc_v0(bp_id) * known_atoms_per_crystal(bp_id))
 
! ELASTIC CONSTANTS

! FIT POLY
rn = 1
DO ec_n = 1, 9
  rn = rn + 1
  a = bp_cn(bp_id, rn, 1)
  b = bp_cn(bp_id, rn, 2)
  s = b - a + 1  
  calc_energies(bp_id, ec_n+1, 1:s) = config_energy(a:b, 6)
  calc_sizes(bp_id, ec_n+1) = s
  CALL fit(bp_data_r(a:b, 3), config_energy(a:b, 6), 2, calc_fitting(bp_id, ec_n, 1:3))
  DO n = 1, s
    CALL poly_calc(calc_fitting(bp_id, ec_n, 1:3), calc_strains(bp_id, ec_n+1, n), calc_energies_fit(bp_id, ec_n+1, n))
  END DO  
END DO


calc_ec(bp_id,1:6,1:6) = 0.0D0
c11 = (2.0D0 / calc_v0(bp_id)) * calc_fitting(bp_id, 1, 3)
c22 = (2.0D0 / calc_v0(bp_id)) * calc_fitting(bp_id, 2, 3)
c33 = (2.0D0 / calc_v0(bp_id)) * calc_fitting(bp_id, 3, 3)
c44 = (1.0D0 / (2.0D0 * calc_v0(bp_id))) * calc_fitting(bp_id, 4, 3)
c55 = (1.0D0 / (2.0D0 * calc_v0(bp_id))) * calc_fitting(bp_id, 5, 3)
c66 = (1.0D0 / (2.0D0 * calc_v0(bp_id))) * calc_fitting(bp_id, 6, 3)
c12 = ((c11 + c22) / 2.0D0) - (calc_fitting(bp_id, 7, 3) / calc_v0(bp_id))
c13 = ((c11 + c33) / 2.0D0) - (calc_fitting(bp_id, 8, 3) / calc_v0(bp_id))
c23 = ((c22 + c33) / 2.0D0) - (calc_fitting(bp_id, 9, 3) / calc_v0(bp_id))

calc_ec(bp_id, 1, 1) = C11
calc_ec(bp_id, 2, 2) = C22
calc_ec(bp_id, 3, 3) = C33
calc_ec(bp_id, 4, 4) = C44
calc_ec(bp_id, 5, 5) = C55
calc_ec(bp_id, 6, 6) = C66

calc_ec(bp_id, 1, 2) = C12
calc_ec(bp_id, 1, 3) = C13
calc_ec(bp_id, 2, 3) = C23

calc_ec(bp_id, 2, 1) = C12
calc_ec(bp_id, 3, 1) = C13
calc_ec(bp_id, 3, 2) = C23

calc_ec_gpa(bp_id,1:6,1:6) = 160.230732254D0 * calc_ec(bp_id,1:6,1:6)

CALL inv(calc_ec(bp_id,1:6,1:6), calc_sc(bp_id,1:6,1:6))
CALL inv(calc_ec_gpa(bp_id,1:6,1:6), calc_sc_gpa(bp_id,1:6,1:6))


!# BULK MODULUS
calc_b0_r(bp_id) = 1.0D0 / (&
                   calc_sc(bp_id, 1, 1) + calc_sc(bp_id, 2, 2) + calc_sc(bp_id, 3, 3) + &
                   2.0D0 * (calc_sc(bp_id, 1, 2) + calc_sc(bp_id, 1, 3) + calc_sc(bp_id, 2, 3)))
calc_b0_v(bp_id) = (1.0D0 / 9.0D0) * (&
                   calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 2, 2) + calc_ec(bp_id, 3, 3)) + &
                   (2.0D0 / 9.0D0) * (&
                   calc_ec(bp_id, 1, 2) + calc_ec(bp_id, 1, 3) + calc_ec(bp_id, 2, 3))
calc_b0_avg(bp_id) = 0.5D0 * (calc_b0_r(bp_id) + calc_b0_v(bp_id))

!# SHEAR MODULUS
calc_g_r(bp_id) = 15.0D0 / (&
                    4.0D0 * (calc_sc(bp_id, 1, 1) + calc_sc(bp_id, 2, 2) + calc_sc(bp_id, 3, 3)) &
                  - 4.0D0 * (calc_sc(bp_id, 1, 2) + calc_sc(bp_id, 1, 3) + calc_sc(bp_id, 2, 3)) &
                  + 3.0D0 * (calc_sc(bp_id, 4, 4) + calc_sc(bp_id, 5, 5) + calc_sc(bp_id, 6, 6)) &
                  )
calc_g_v(bp_id) = (1.0D0 / 15.0D0) * (&
                  calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 2, 2) + calc_ec(bp_id, 3, 3) &
                  - calc_ec(bp_id, 1, 2) - calc_ec(bp_id, 1, 3) - calc_ec(bp_id, 2, 3) &
                  ) &
                  + (1.0D0/5.0D0) * (&
                  calc_ec(bp_id, 4, 4) + calc_ec(bp_id, 5, 5) + calc_ec(bp_id, 6, 6) &
                  )
calc_g_avg(bp_id) = 0.5D0 * (calc_g_r(bp_id) + calc_g_v(bp_id))


calc_g_vec(bp_id, 1) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 4, 4))
calc_g_vec(bp_id, 2) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 5, 5))
calc_g_vec(bp_id, 3) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 6, 6))

calc_e(bp_id) = (9 * calc_b0_avg(bp_id) * calc_g_avg(bp_id)) / (3.0D0 * calc_b0_avg(bp_id) + calc_g_avg(bp_id))

calc_e_vec(bp_id, 1) = 1.0D0 / calc_sc(bp_id, 1, 1)
calc_e_vec(bp_id, 2) = 1.0D0 / calc_sc(bp_id, 2, 2)
calc_e_vec(bp_id, 3) = 1.0D0 / calc_sc(bp_id, 3, 3)


calc_v(bp_id) = (3.0D0 * calc_b0_avg(bp_id) * calc_g_avg(bp_id)) / ( 3 * calc_b0_avg(bp_id) + calc_g_avg(bp_id))


calc_v_vec(bp_id, 1, 1) = 1.0D0
calc_v_vec(bp_id, 2, 1) = -(calc_sc(bp_id, 1, 2) * calc_e_vec(bp_id, 2))
calc_v_vec(bp_id, 3, 1) = -(calc_sc(bp_id, 1, 3) * calc_e_vec(bp_id, 3))
calc_v_vec(bp_id, 1, 2) = -(calc_sc(bp_id, 2, 1) * calc_e_vec(bp_id, 1))
calc_v_vec(bp_id, 2, 2) = 1.0D0
calc_v_vec(bp_id, 3, 2) = -(calc_sc(bp_id, 2, 3) * calc_e_vec(bp_id, 3))
calc_v_vec(bp_id, 1, 3) = -(calc_sc(bp_id, 3, 1) * calc_e_vec(bp_id, 1))
calc_v_vec(bp_id, 2, 3) = -(calc_sc(bp_id, 3, 2) * calc_e_vec(bp_id, 2))
calc_v_vec(bp_id, 3, 3) = 1.0D0

b0 = 160.230732254D0 * calc_b0_avg(bp_id)
g = 160.230732254D0 * calc_g_avg(bp_id)


calc_vl(bp_id) = sqrt((b0 + (4.0D0/3.0D0) * g) / calc_rho(bp_id))
calc_vt(bp_id) = sqrt(g / calc_rho(bp_id))
calc_vm(bp_id) = ((1.0D0/3.0D0) * (2.0D0 / (calc_vt(bp_id)**3)) + 1.0D0 / (calc_vl(bp_id)**3))**(-(1.0D0/3.0D0))

calc_melting(bp_id) =   598.0D0 &
                      + 6.66D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3)) &
                      - 0.003D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3))**2  
       
h = 6.62607004D-34
k = 1.38064852D-23
pi = 3.1415926535898D0
na = 6.02214076D23
calc_debye(bp_id) = ((h / k) * &
                    ((( 3.0D0 * known_atoms_per_crystal(bp_id)) / (4.0D0 * pi)) * &
                    ((na * calc_rho(bp_id)) / (known_amu_per_crystal(bp_id) * 1.66D-27))**(1.0D0/3.0D0)) &
                    * calc_vm(bp_id))


! FINALLY CALCULATE RSS
rss = 0.0D0
rss_total = 0.0D0
rss_w = 0.0D0
rss_total_w = 0.0D0

IF(known_set(bp_id, 1))THEN
  rss(bp_id, 1) = (known_alat(bp_id) - calc_alat(bp_id))**2
END IF
IF(known_set(bp_id, 2))THEN
  rss(bp_id, 2) = (known_e0(bp_id) - calc_e0(bp_id))**2
END IF
IF(known_set(bp_id, 3))THEN
  rss(bp_id, 3) = (known_b0(bp_id) - calc_b0(bp_id))**2
END IF
IF(known_set(bp_id, 4))THEN
  rss(bp_id, 4) = SUM((known_ec(bp_id,1:6,1:6) - calc_ec(bp_id,1:6,1:6))**2)
END IF
IF(known_set(bp_id, 5))THEN
  rss(bp_id, 5) = (known_g(bp_id) - calc_g_avg(bp_id))**2
END IF
IF(known_set(bp_id, 6))THEN
  rss(bp_id, 6) = (known_e(bp_id) - calc_e(bp_id))**2
END IF
IF(known_set(bp_id, 7))THEN
  rss(bp_id, 7) = (known_v(bp_id) - calc_v(bp_id))**2
END IF
rss_total(bp_id) = SUM(rss(bp_id, 1:7))
rss_total_rss = rss_total_rss + rss_total(bp_id)

rss_w(bp_id, 1:7) = rss_weights(1:7) * rss(bp_id, 1:7)
rss_total_w(bp_id) = SUM(rss_w(bp_id, 1:7))
rss_total_rss_w = rss_total_rss_w + rss_total_w(bp_id)

END SUBROUTINE calculate_bp_inner