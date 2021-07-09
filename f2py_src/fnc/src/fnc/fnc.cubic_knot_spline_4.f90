! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
!
!  Fixed parameters:  x positions, zbl a, zbl b, zbl mult, end point x, y(x), y'(x)
!  Variable parameters: y(x) y'(x)
!
!    test.fn = 'cubic_knot_spline_4'
!    pf = 0.0, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.5, 4.0,  4.2,  5.0, 6.5, 0.0, 0.0, 26.0, 26.0, 1.0
!    p = 0.0, 1.0, 1.1, 1.0, 0.6, 0.5, 0.2, 0.23, 4.7, 0.21, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0


SUBROUTINE cubic_knot_spline_4(r, p, p_fixed, y)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: n
REAL(kind=DoubleReal) :: nodes(1:size(p) / 2 + 1,1:4)
INTEGER(kind=StandardInteger) :: node_count
REAL(kind=DoubleReal) :: qa, qb, zbl_mult, grad_soften
INTEGER(kind=StandardInteger) :: psize
LOGICAL :: complete
!############################################################

! P and PF must be the same length
! PF contains x value of node
! P contains y value of node
nodes = 0.0D0

psize = size(p,1) / 2
node_count = size(p,1) / 2 + 1

nodes(1:psize, 1) = p_fixed(1:psize)
nodes(1:psize, 2) = p(1:psize)
nodes(1:psize, 3) = p(psize+1:psize+psize)

qa = p_fixed(psize + 4)
qb = p_fixed(psize + 5)
zbl_mult = p_fixed(psize + 6)

nodes(node_count, 1) = p_fixed(psize + 1)
nodes(node_count, 2) = p_fixed(psize + 2)
nodes(node_count, 3) = p_fixed(psize + 3)

IF(r .LT. nodes(1, 1))THEN   
  y = nodes(1, 2)
ELSE IF(r .GT. nodes(node_count, 1))THEN  
  y = nodes(node_count, 2)
ELSE  
  ! Add ZBL
  y = y + 1.0D0 * f_zbl_ackland_mendelev(r, qa, qb) 
  
  ! IF EXACTLY A NODE
	complete = .FALSE.
	DO n = 1, node_count  
	 	IF(r .EQ. nodes(n,1))THEN      
	  	complete = .TRUE.
	  	y = y + nodes(n,2)
	 	END IF
	END DO
  
  IF(complete .EQV. .FALSE.)THEN
    y = y + f_cubic_knot_spline(r, nodes)
  END IF

END IF

END SUBROUTINE cubic_knot_spline_4


! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_4_v(r, p, p_fixed, y)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
! Loop through all the values in r(:), calculate and store in y(:)
DO n = 1, SIZE(r,1)
  CALL cubic_knot_spline_4(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_4_v















