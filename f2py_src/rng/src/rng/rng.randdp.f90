SUBROUTINE randdp(ra, rb, rn)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: ra
REAL(kind=DoubleReal), INTENT(IN) :: rb
REAL(kind=DoubleReal), INTENT(OUT) :: rn
!#############################################################
REAL(kind=DoubleReal) :: rdp
!#############################################################
CALL RANDOM_NUMBER(rdp) 
rn = ra + (rb - ra) * rdp
!#############################################################
END SUBROUTINE randdp