SUBROUTINE randint_2d(ra, rb, a, b, rn)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: ra
INTEGER(KIND=StandardInteger), INTENT(IN) :: rb
INTEGER(KIND=StandardInteger), INTENT(IN) :: a
INTEGER(KIND=StandardInteger), INTENT(IN) :: b
INTEGER(KIND=StandardInteger), INTENT(OUT) :: rn(1:a, 1:b)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
REAL(kind=DoubleReal) :: rdp(1:a, 1:b)
!#############################################################
CALL RANDOM_NUMBER(rdp(:,:)) 
rn(:,:) = MIN(ra + FLOOR((rb - ra + 1) * rdp(:,:)), rb)
!#############################################################
END SUBROUTINE randint_2d
