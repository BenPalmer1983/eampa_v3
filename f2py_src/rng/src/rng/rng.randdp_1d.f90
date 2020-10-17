SUBROUTINE randdp_1d(ra, rb, a, rn)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: ra
REAL(kind=DoubleReal), INTENT(IN) :: rb
INTEGER(KIND=StandardInteger), INTENT(IN) :: a
REAL(kind=DoubleReal), INTENT(OUT) :: rn(1:a)
!#############################################################
REAL(kind=DoubleReal) :: rdp(1:a)
!#############################################################
CALL RANDOM_NUMBER(rdp) 
rn(:) = ra + (rb - ra) * rdp(:)
!#############################################################
END SUBROUTINE randdp_1d