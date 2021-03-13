SUBROUTINE randdp_1d(ra, rb, a, b, rn)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: ra
REAL(kind=DoubleReal), INTENT(IN) :: rb
INTEGER(KIND=StandardInteger), INTENT(IN) :: a
INTEGER(KIND=StandardInteger), INTENT(IN) :: b
REAL(kind=DoubleReal), INTENT(OUT) :: rn(1:a, 1:b)
!#############################################################
REAL(kind=DoubleReal) :: rdp(1:a, 1:b)
!#############################################################
CALL RANDOM_NUMBER(rdp) 
rn(:,:) = ra + (rb - ra) * rdp(:,:)
!#############################################################
END SUBROUTINE randdp_1d