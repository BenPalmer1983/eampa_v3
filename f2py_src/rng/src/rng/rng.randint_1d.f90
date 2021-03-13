SUBROUTINE randint_1d(ra, rb, a, rn)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: ra
INTEGER(KIND=StandardInteger), INTENT(IN) :: rb
INTEGER(KIND=StandardInteger), INTENT(IN) :: a
INTEGER(KIND=StandardInteger), INTENT(OUT) :: rn(1:a)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
REAL(kind=DoubleReal) :: rdp(1:a)
!#############################################################
CALL RANDOM_NUMBER(rdp(:)) 
rn(:) = MIN(ra + FLOOR((rb - ra + 1) * rdp(:)), rb)
!#############################################################
END SUBROUTINE randint_1d
