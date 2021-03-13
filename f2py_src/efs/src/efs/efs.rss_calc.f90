
SUBROUTINE rss_calc()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn, n
INTEGER(KIND=StandardInteger) :: a, b, rn
REAL(kind=DoubleReal) :: t, f
REAL(kind=DoubleReal) :: we, wf, ws
REAL(kind=DoubleReal) :: energy_rss
REAL(kind=DoubleReal) :: force_rss(1:3)
REAL(kind=DoubleReal) :: stress_rss(1:9)

!###########################################################

CALL start_t()

max_density = 0.0D0

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




residuals_size = 0
DO cn = 1, cc
  residuals_size = residuals_size + 1   ! energy
  IF(key(cn, 12) .EQ. 1)THEN
    residuals_size = residuals_size + 3   ! force
  END IF
  IF(key(cn, 13) .EQ. 1)THEN
    residuals_size = residuals_size + 9   ! stress
  END IF
END DO

! MAKE RESIDUALS ARRAY - WEIGHTED RSS ONLY
IF(ALLOCATED(residuals))THEN
  DEALLOCATE(residuals)
END IF
ALLOCATE(residuals(1:residuals_size))
residuals(1:residuals_size) = 0.0D0
config_rss = 0.0D0


we = rss_weights(1) * rss_weights(2)
wf = rss_weights(1) * rss_weights(3)
ws = rss_weights(1) * rss_weights(4)

rn = 1
DO cn = 1, cc
  ! Energy
  energy_rss = (energies(cn) - config_energy(cn, 3))**2
  config_rss(cn, 1) = energy_rss
  config_rss(cn, 5) = we * energy_rss
  residuals(rn) = energy_rss
  rn = rn + 1

  IF(key(cn, 12) .EQ. 1)THEN
    !config_rss(cn, 2) = SUM((config_forces(key(cn, 1):key(cn, 2), 1:3) - forces(key(cn, 1):key(cn, 2), 1:3))**2)
    force_rss(1) = SUM((config_forces(cn, 1:key(cn, 20), 1) - forces(cn, 1:key(cn, 20), 1)))
    force_rss(2) = SUM((config_forces(cn, 1:key(cn, 20), 1) - forces(cn, 1:key(cn, 20), 1)))
    force_rss(3) = SUM((config_forces(cn, 1:key(cn, 20), 1) - forces(cn, 1:key(cn, 20), 1)))

    config_rss(cn, 2) = SUM(force_rss**2)
    config_rss(cn, 6) = wf * SUM(force_rss**2)
    
    residuals(rn:rn+2) = force_rss(1:3)

    rn = rn + 3
  END IF


  IF(key(cn, 13) .EQ. 1)THEN
    stress_rss(1:3) = (config_stresses(cn, 1, :) - stresses(cn, 1, :))
    stress_rss(4:6) = (config_stresses(cn, 2, :) - stresses(cn, 2, :))
    stress_rss(7:9) = (config_stresses(cn, 3, :) - stresses(cn, 3, :))

    config_rss(cn, 3) = SUM(stress_rss**2)
    config_rss(cn, 7) = ws * SUM(stress_rss**2)

    residuals(rn:rn+8) = stress_rss(1:9)

    rn = rn + 9
  END IF

  config_rss(cn, 4) = SUM(config_rss(cn, 1:3))
  config_rss(cn, 8) = SUM(config_rss(cn, 5:7))

END DO




! TOTALS OVER ALL CONFIGS
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













