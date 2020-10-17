SUBROUTINE calculate_es()
!###########################################################
INTEGER(KIND=StandardInteger) :: calc_id
INTEGER(KIND=StandardInteger) :: cn, n
!###########################################################

! Reset total rss over all 
rss_total_rss = 0.0D0
rss_total_rss_w = 0.0D0

! CALC ENERGIES
CALL energy()

! Loop through
calc_id = 1
DO WHILE(calc_keys_i(calc_id, 1) .NE. -1)
  cn = calc_keys_i(calc_id, 3)
  surface_energy(calc_id) = (config_energy(cn+1, 3) - config_energy(cn, 3)) / (2.0D0 * surface_area(calc_id, 1))
  surface_energy_mjmsq(calc_id) = 1602.0D0 * surface_energy(calc_id)
  !print *, calc_id, surface_energy(calc_id), surface_energy_mjmsq(calc_id)
  calc_id = calc_id + 1
END DO



!CALL force()
!print *, key(1, 1), key(1, 2)
!DO n=key(1, 1), key(1, 2)
!  print *, config_forces(n, 1, :)
!END DO

END SUBROUTINE calculate_es
