

SUBROUTINE search_f(f_type, option_a, option_b, option_g, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b, option_g
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
LOGICAL :: in_cache
!###########################################################
! Default Out
fx = 0.0D0
IF(function_cache)THEN
  ! Get Key   
  CALL get_pot_key(f_type, option_a, option_b, option_g, pkey)
  IF(pkey .GT. 0)THEN
    IF(pot_f_type(pkey) .GT. 0)THEN
      CALL get(pkey, x, in_cache, fx)
      IF(.NOT. in_cache)THEN
        CALL f(pot_f_name(pkey), x, params(pkey, 1:p_count(pkey)), &
            paramsfixed(pkey, 1:pf_count(pkey)), fx)
        CALL set(pkey, x, fx)
      END IF
    END IF
  END IF
ELSE
  CALL get_pot_key(f_type, option_a, option_b, option_g, pkey)
  CALL f(pot_f_name(pkey), x, params(pkey, 1:p_count(pkey)), &
    paramsfixed(pkey, 1:pf_count(pkey)), fx)
END IF
END SUBROUTINE search_f





SUBROUTINE search_grad(f_type, option_a, option_b, option_g, x, fx)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: f_type, option_a, option_b, option_g
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: fx
!###########################################################
INTEGER(kind=StandardInteger) :: pkey
LOGICAL :: in_cache
!###########################################################
! Default Out
fx = 0.0D0
! Get Key
CALL get_pot_key(f_type, option_a, option_b, option_g, pkey)
IF(pkey .GT. 0)THEN
  IF(pot_f_type(pkey) .GT. 0)THEN
    CALL get(pkey + key_max, x, in_cache, fx)
    IF(.NOT. in_cache)THEN
      CALL fgrad(pot_f_name(pkey), x, params(pkey, 1:p_count(pkey)), &
             paramsfixed(pkey, 1:p_count(pkey)), fx)
      CALL set(pkey + key_max, x, fx)
    END IF
  END IF
END IF
END SUBROUTINE search_grad









