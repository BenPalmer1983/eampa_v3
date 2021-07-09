



! KEYS
SUBROUTINE get_pot_key(ftype, a, b, group, key)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: ftype, a, b, group
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
!###########################################################
IF(ftype .EQ. 1)THEN
  CALL pair_key(a, b, key)
ELSE IF(ftype .EQ. 2)THEN
  CALL dens_key(a, b, group, key)
ELSE IF(ftype .EQ. 3)THEN
  CALL embe_key(a, group, key)
END IF

END SUBROUTINE get_pot_key





! KEYS
SUBROUTINE pair_key(a, b, key)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, b
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
! Min/Max
key = 0
IF(a .GT. label_max .OR. b .GT. label_max)THEN
  key = 0
ELSE
  min_key = a
  max_key = b
  IF(a .GT. b)THEN
    max_key = a
    min_key = b
  END IF
  key = (max_key*(max_key-1))/2 + min_key
END IF
END SUBROUTINE pair_key




! KEYS
SUBROUTINE dens_key(a, b, fgroup, key)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, b, fgroup
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
key = 0
IF(a .GT. label_max .OR. a .GT. label_max .OR. fgroup .GT. fgroup_max)THEN
  key = 0
ELSE
  IF(b .EQ. 0)THEN
    key = dens_key_offset + (fgroup - 1) * (label_max + pair_max) + a
  ELSE IF(b .GT. 0)THEN
    min_key = a
    max_key = b
    IF(a .GT. b)THEN
      max_key = a
      min_key = b
    END IF
    key = (max_key*(max_key-1))/2 + min_key
    key = dens_key_offset + (fgroup - 1) * (label_max + pair_max) + label_max + &
          (max_key*(max_key-1))/2 + min_key
  END IF
END IF
END SUBROUTINE dens_key


! KEYS
SUBROUTINE embe_key(a, fgroup, key)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, fgroup
INTEGER(kind=StandardInteger), INTENT(OUT) :: key
INTEGER(kind=StandardInteger) :: min_key, max_key
!###########################################################
key = 0
IF(a .GT. label_max .OR. a .GT. label_max .OR. fgroup .GT. fgroup_max)THEN
  key = 0
ELSE
  key = embe_key_offset + (fgroup - 1) * (label_max + pair_max) + a
END IF
END SUBROUTINE embe_key









