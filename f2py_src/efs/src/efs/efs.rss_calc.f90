
SUBROUTINE rss_calc()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn, n
INTEGER(KIND=StandardInteger) :: a, b
REAL(kind=DoubleReal) :: t
!###########################################################

CALL start_t()


! CALCULATE E, F, S
!$OMP PARALLEL
!$OMP DO
DO cn = 1, cc
  IF(key(cn, 11) .EQ. 1 .AND. key(cn, 12) .EQ. 0)THEN
    CALL energy_inner(cn)    
  ELSE IF(key(cn, 11) .EQ. 1 .AND. key(cn, 12) .EQ. 1 .AND. key(cn, 13) .EQ. 0)THEN
    CALL energy_force_inner(cn)
  ELSE IF((key(cn, 11) .EQ. 1 .AND. key(cn, 13) .EQ. 1))THEN
    CALL energy_force_stress_inner(cn)
  END IF
END DO  
!$OMP END DO
!$OMP END PARALLEL


! CALCULATE RSS
config_rss = 0.0D0
DO cn = 1, cc
  config_rss(cn, 1) = (energies(cn) - config_energy(cn, 3))**2
  IF(key(cn, 12) .EQ. 1)THEN
    config_rss(cn, 2) = SUM((config_forces(key(cn, 1):key(cn, 2), 1:3) - forces(key(cn, 1):key(cn, 2), 1:3))**2)
  END IF
  IF(key(cn, 13) .EQ. 1)THEN
    config_rss(cn, 3) = SUM((config_stresses(3 * (cn-1) + 1: 3 * (cn-1) + 3,:) - stresses(3 * (cn-1) + 1: 3 * (cn-1) + 3,:))**2)
  END IF
  config_rss(cn, 4) = SUM(config_rss(cn, 1:3)) 
  
  config_rss(cn, 5) = rss_weights(1) * rss_weights(2) * config_rss(cn, 1)   ! ENERGY
  config_rss(cn, 6) = rss_weights(1) * rss_weights(3) * config_rss(cn, 2)   ! FORCE
  config_rss(cn, 7) = rss_weights(1) * rss_weights(4) * config_rss(cn, 3)   ! STRESS
  config_rss(cn, 8) = SUM(config_rss(cn, 5:7)) 
  
  !print *, config_rss(cn, 1), config_energy(cn, 3), energies(cn)
  !print *, rss_weights(2), rss_weights(3), rss_weights(4)
END DO

energy_rss = SUM(config_rss(1:cc, 1))
force_rss = SUM(config_rss(1:cc, 2))
stress_rss = SUM(config_rss(1:cc, 3))
total_rss = SUM(config_rss(1:cc, 4))
energy_rss_weighted = SUM(config_rss(1:cc, 5))
force_rss_weighted = SUM(config_rss(1:cc, 6))
stress_rss_weighted = SUM(config_rss(1:cc, 7))
total_rss_weighted = SUM(config_rss(1:cc, 8))





CALL end_t(rss_timer)
rss_timer_sum = rss_timer_sum + rss_timer
END SUBROUTINE rss_calc













