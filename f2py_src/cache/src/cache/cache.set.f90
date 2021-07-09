SUBROUTINE set(f_key, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_key
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(IN) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: i, n, m
LOGICAL :: in_cache
!###########################################################
IF(.NOT. isnan(x) .AND. (abs(x) .LT. 1.0D300))THEN
  IF(x .EQ. 0.0D0)THEN
    n = 1 
  ELSE
    n = mod(floor(a_size * (abs(x) / (10**(1.0D0*floor(log10(abs(x))))))), a_size) + 1
  END IF
  in_cache = .FALSE.
  m = 0
  DO i = 1, cache_count(f_key, n)
    IF(cache_data(f_key, n, i, 1) .EQ. x)THEN
      in_cache = .TRUE.
      m = i
      EXIT
    END IF
  END DO

  IF((in_cache .EQV. .FALSE.) .AND. (cache_count(f_key, n) .LT. b_size))THEN
    cache_count(f_key, n) = cache_count(f_key, n) + 1
    cache_data(f_key, n, cache_count(f_key, n), 1) = x
    cache_data(f_key, n, cache_count(f_key, n), 2) = fx
  END IF
END IF
END SUBROUTINE set