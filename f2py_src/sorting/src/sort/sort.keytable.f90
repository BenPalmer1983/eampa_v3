SUBROUTINE make_keytable(n)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
INTEGER(KIND=StandardInteger) :: i
!#############################################################

IF(ALLOCATED(keytable))THEN
  DEALLOCATE(keytable)
END IF
ALLOCATE(keytable(1:n))

DO i = 1, n
  keytable(i) = i
END DO

END SUBROUTINE make_keytable

