SUBROUTINE init()
!###########################################################
IF(initialised .EQV. .FALSE.)THEN
  initialised = .TRUE.

  IF(ALLOCATED(cache_count)) DEALLOCATE(cache_count)
  IF(ALLOCATED(cache_data)) DEALLOCATE(cache_data)
  ALLOCATE(cache_count(1:f_keys, 1:a_size))
  ALLOCATE(cache_data(1:f_keys, 1:a_size, 1:b_size, 1:2))

  cache_count = 0
  cache_data = 0.0D0
END IF

END SUBROUTINE init