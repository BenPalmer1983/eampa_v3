! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
!
! PF   Qa  Qb  zbl_mult  rs   v(rs)  v'(rs)  re  v(re)  v'(re)
! P    node r    node v(r)


SUBROUTINE cubic_knot_spline_3(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb, yzbl
REAL(kind=DoubleReal) :: start_x, start_y, start_yp
REAL(kind=DoubleReal) :: end_x, end_y, end_yp, zbl_mult, grad_soften
INTEGER(kind=StandardInteger) :: i, j, k, n, n_p, n_pf
REAL(kind=DoubleReal) :: nodes(1:size(p) / 2 + 2,1:4)
INTEGER(kind=StandardInteger) :: node_count
REAL(kind=DoubleReal) :: qa, qb, rzbl
REAL(kind=DoubleReal) :: start(1:4) = 0.0D0
REAL(kind=DoubleReal) :: end(1:4) = 0.0D0
INTEGER(kind=StandardInteger) :: psize
INTEGER(kind=StandardInteger) :: na, nb, nsize
INTEGER(kind=StandardInteger) :: interp_size = 5
LOGICAL :: complete
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
!############################################################

! P and PF must be the same length
! PF contains x value of node
! P contains y value of node
nodes = 0.0D0

psize = size(p,1) / 2

node_count = size(p) / 2 + 2
nodes(2:psize+1, 1) = p(1:psize)
nodes(2:psize+1, 2) = p(psize+1:psize+psize)

qa = p_fixed(1)
qb = p_fixed(2)
zbl_mult = p_fixed(3)

nodes(1, 1) = p_fixed(4)
nodes(1, 2) = p_fixed(5)
nodes(1, 3) = p_fixed(6)  ! Re Set Later

nodes(node_count, 1) = p_fixed(7)
nodes(node_count, 2) = p_fixed(8)
nodes(node_count, 3) = p_fixed(9)  ! Re Set Later

grad_soften = p_fixed(10)

IF(r .LT. nodes(1, 1))THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = 0.0D0
ELSE IF(r .GT. nodes(node_count, 1))THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = 0.0D0
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
    CALL interp_fill_derivatives(nodes(1:node_count, 1:4), 4)
    nodes(1, 3) = p_fixed(6)
    nodes(1:node_count, 3) = grad_soften * nodes(1:node_count, 3)
    nodes(node_count, 3) = p_fixed(9)
    y = y + f_cubic_knot_spline(r, nodes)
  END IF

END IF

END SUBROUTINE cubic_knot_spline_3


! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_3_v(r, p, p_fixed, y)
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
  CALL cubic_knot_spline_3(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_3_v















