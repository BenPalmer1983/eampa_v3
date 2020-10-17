SUBROUTINE randint(ra, rb, rn)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: ra
INTEGER(KIND=StandardInteger), INTENT(IN) :: rb
INTEGER(KIND=StandardInteger), INTENT(OUT) :: rn
!#############################################################
REAL(kind=DoubleReal) :: rdp
!#############################################################
CALL RANDOM_NUMBER(rdp) 
rn = MIN(ra + FLOOR((rb - ra + 1) * rdp), rb)
!#############################################################
END SUBROUTINE randint
