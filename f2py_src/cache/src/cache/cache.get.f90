SUBROUTINE get(f_key, x, in_cache, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_key
REAL(kind=DoubleReal), INTENT(IN) :: x
LOGICAL, INTENT(OUT) :: in_cache
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: i, n
!###########################################################
CALL init()
in_cache = .FALSE.
IF(.NOT. isnan(x))THEN
fx = 0.0D0
IF(x .EQ. 0.0D0)THEN
  n = 1
ELSE
  n = mod(floor(a_size * (abs(x) / (10**(1.0D0*floor(log10(abs(x))))))), a_size) + 1
END IF
DO i = 1, cache_count(f_key, n)
  IF(cache_data(f_key, n, i, 1) .EQ. x)THEN
    in_cache = .TRUE.
    fx = cache_data(f_key, n, i, 2)
    EXIT
  END IF
END DO
END IF
END SUBROUTINE get