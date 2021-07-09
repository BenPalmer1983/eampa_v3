

SUBROUTINE fzbl(r, qa, qb, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: qa
REAL(kind=DoubleReal), INTENT(IN) :: qb
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: rs
REAL(kind=DoubleReal) :: ep
!############################################################
IF(r .LE. 0.0D0)THEN
  y = 1.0D10
ELSE
  rs = 0.4683766 / (qa**(2.0D0/3.0D0)+qb**(2.0D0/3.0D0))
  ep = 0.1818D0 * exp(-3.2D0 * (r/rs)) 
  ep = ep + 0.5099D0 * exp(-0.9423D0 * (r/rs)) 
  ep = ep + 0.2802 * exp(-0.4029D0 * (r/rs)) 
  ep = ep + 0.02817 * exp(-0.2016D0 * (r/rs)) 
  y = ((qa * qb) / r) * ep
END IF
END SUBROUTINE fzbl




SUBROUTINE fzbl_v(r, qa, qb, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: qa
REAL(kind=DoubleReal), INTENT(IN) :: qb
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL fzbl(r(n), qa, qb, y(n))
END DO
END SUBROUTINE fzbl_v







SUBROUTINE fzbl_dydr(r, qa, qb, dydr)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: qa
REAL(kind=DoubleReal), INTENT(IN) :: qb
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: rs
REAL(kind=DoubleReal) :: ep
REAL(kind=DoubleReal) :: epp
!############################################################
IF(r .LE. 0.0D0)THEN
  dydr = -1.0D10
ELSE
  rs = 0.4683766 / (qa**(2.0D0/3.0D0)+qb**(2.0D0/3.0D0))
  
  
  ep =      0.1818D0 * exp(r * (-3.2D0/rs)) 
  ep = ep + 0.5099D0 * exp(r * (-0.9423D0/rs)) 
  ep = ep + 0.2802 * exp(r * (-0.4029D0/rs)) 
  ep = ep + 0.02817 * exp(r * (-0.2016D0/rs)) 
  
  
  epp =      (-3.2D0/rs) * 0.1818D0 * exp(r * (-3.2D0/rs)) 
  epp = ep + (-0.9423D0/rs) * 0.5099D0 * exp(r * (-0.9423D0/rs)) 
  epp = ep + (-0.4029D0/rs) * 0.2802 * exp(r * (-0.4029D0/rs)) 
  epp = ep + (-0.2016D0/rs) * 0.02817 * exp(r * (-0.2016D0/rs)) 
  
  dydr = ((qa * qb) / r) * epp - ((qa * qb) / (r**2)) * ep
END IF
END SUBROUTINE fzbl_dydr



SUBROUTINE fzbl_dydr_v(r, qa, qb, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: qa
REAL(kind=DoubleReal), INTENT(IN) :: qb
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL fzbl_dydr(r(n), qa, qb,  y(n))
END DO
END SUBROUTINE fzbl_dydr_v

