! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
!
!  Fixed parameters:  x positions, zbl a, zbl b, zbl mult, end point x, y(x), y'(x)
!  Variable parameters: y(x) y'(x)
!
!    test.fn = 'cubic_knot_spline_4'
!         zbla  zblb  zblm xstart  xend  v(xend)  v'(xend)
!    pf = 26.0, 26.0, 1.0, 0.0,    6.5,  0.0,     0.0 
!    p =  1.0  2.0  3.0  4.0  5.0      0.2  0.3  -0.1  0.01  0.001   0.0  0.0  0.0  0.0  0.0   0.0   0.0    (last two v(xstart) v'(xstart))


SUBROUTINE cubic_knot_spline_5(r, p, p_fixed, y)
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
INTEGER(kind=StandardInteger) ::n
REAL(kind=DoubleReal) :: nodes(1:(size(p,1)-2) / 3 + 2,1:4)
REAL(kind=DoubleReal) :: qa, qb, zbl_mult
INTEGER(kind=StandardInteger) :: psize, pfsize, node_count, pstride
LOGICAL :: complete
!############################################################
psize = SIZE(p, 1)
pfsize = SIZE(p_fixed, 1)
pstride = (size(p,1)-2) / 3
node_count = pstride + 2
nodes(:,:) = 0.0D0
nodes(2:node_count-1, 1) = p(1:pstride)
nodes(2:node_count-1, 2) = p(pstride+1:2*pstride)
nodes(2:node_count-1, 3) = p(2*pstride+1:3*pstride)

qa = p_fixed(1)
qb = p_fixed(2)
zbl_mult = p_fixed(3)
nodes(1, 1) = p_fixed(4)
nodes(1, 2) = p(3*pstride+1)
nodes(1, 3) = p(3*pstride+2)
nodes(node_count, 1) = p_fixed(5)
nodes(node_count, 2) =  p_fixed(6)
nodes(node_count, 3) =  p_fixed(7)
!print *, pstride, node_count
!print *, nodes(:, 1)
!print *, nodes(:, 2)
!print *, nodes(:, 3)
!print *, ""
IF(r .LT. nodes(1, 1))THEN   
  y =  nodes(1, 2)
ELSE IF(r .GT. nodes(node_count, 1))THEN  
  y =  nodes(node_count, 2)
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


END SUBROUTINE cubic_knot_spline_5


! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_5_v(r, p, p_fixed, y)
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
  CALL cubic_knot_spline_5(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_5_v















