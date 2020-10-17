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
  !CALL calculate_bp_inner(bp_id)
  print *, calc_id
  print *, calc_keys_i(calc_id, 3), calc_keys_i(calc_id, 4)
  DO cn = calc_keys_i(calc_id, 3), calc_keys_i(calc_id, 4)
    print *,"   ",cn, config_energy(cn, 3)
  END DO
  print *, calc_keys_r(calc_id,3)  ! Area
  
  calc_id = calc_id + 1
END DO




!CALL force()
!print *, key(1, 1), key(1, 2)
!DO n=key(1, 1), key(1, 2)
!  print *, config_forces(n, 1, :)
!END DO

END SUBROUTINE calculate_es
