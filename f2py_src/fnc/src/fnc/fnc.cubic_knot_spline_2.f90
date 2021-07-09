! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
! Pfixed    q1, q2, zbl_cutoff  node1 node 2...node n end_node



SUBROUTINE cubic_knot_spline_2(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: nodes_in(1:SIZE(p,1)+1,1:4)
REAL(kind=DoubleReal) :: nodes(1:10000,1:4)
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

psize = SIZE(p,1)

nodes_in(:,:) = 0.0D0
nodes_in(1:psize, 1) = p_fixed(1:psize)
nodes_in(1:psize, 2) = p(1:psize)
nodes_in(psize+1, 1) = p_fixed(psize+1)
nodes_in(psize+1, 2) = p_fixed(psize+2)
nodes_in(psize+1, 3) = p_fixed(psize+3)

qa = p_fixed(psize + 4)
qb = p_fixed(psize + 5)
zbl_mult = p_fixed(psize + 6)
grad_soften = p_fixed(psize + 7)

IF(SIZE(p_fixed,1) .LT. SIZE(p,1) + 7)THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = 0.0D0
ELSE IF(r .GT. nodes_in(psize+1, 1))THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = nodes_in(psize+1, 2)
ELSE  
  ! Add ZBL
  y = y + 1.0D0 * f_zbl_ackland_mendelev(r, qa, qb) 
  
  ! IF EXACTLY A NODE
	complete = .FALSE.
	DO n = 1, SIZE(nodes_in, 1)  
	 	IF(r .EQ. nodes_in(n,1))THEN      
	  	complete = .TRUE.
	  	y = y + nodes_in(n,2)
	 	END IF
	END DO
  
  IF(complete .EQV. .FALSE.)THEN
    CALL interp_fill_derivatives(nodes_in(1:psize+1, 1:4), 4)
    nodes_in(1:psize+1, 3) = grad_soften * nodes_in(1:psize+1, 3)
    y = y + f_cubic_knot_spline(r, nodes_in)
  END IF

END IF

END SUBROUTINE cubic_knot_spline_2


! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_2_v(r, p, p_fixed, y)
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
  CALL cubic_knot_spline_2(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_2_v















