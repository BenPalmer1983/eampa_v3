! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
! Pfixed    q1, q2, zbl_cutoff  node1 node 2...node n end_node



SUBROUTINE cubic_knot_spline_fixed_end_pair(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb
REAL(kind=DoubleReal) :: start_x, start_y, start_yp
REAL(kind=DoubleReal) :: end_x, end_y, end_yp
INTEGER(kind=StandardInteger) :: i, j, k, n, n_p, n_pf
REAL(kind=DoubleReal) :: nodes_in(1:SIZE(p,1)+2,1:4)
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

IF(SIZE(p_fixed,1) .LT. SIZE(p,1) + 6)THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = 0.0D0
ELSE IF(r .GT. p_fixed(psize + 4))THEN   !  qa, qb, rzbl, r, v(r), v'(r) ..... nsize
  y = p_fixed(psize + 5)
ELSE
  
  qa = p_fixed(psize + 1)
  qb = p_fixed(psize + 2)
  rzbl = p_fixed(psize + 3)
  end(1) = p_fixed(psize + 4)
  end(2) = p_fixed(psize + 5)
  end(3) = p_fixed(psize + 6)

  IF(r .LT. rzbl)THEN
    CALL fzbl(r, qa, qb, y)
  ELSE
    nsize = 1
    IF((SIZE(p_fixed,1) .EQ. psize + 7))THEN
  	  nsize = INT(CEILING(p_fixed(SIZE(p, 1)+7))) 
    END IF

    nodes_in = 0.0D0
    y = 0.0D0

    start(1) = rzbl
    CALL fzbl(rzbl, qa, qb, start(2))
    CALL fzbl_dydr(rzbl, qa, qb, start(3))


    ! LOAD X and Y POINTS
    nodes_in(1,1) = start(1)
    nodes_in(1,2) = start(2)
    nodes_in(1,3) = start(3)
    nodes_in(1,4) = 0.0D0

    nodes_in(2:psize,1) = p_fixed(1:psize)      ! r
    nodes_in(2:psize,2) = p(:)                  ! v(r)
    nodes_in(2:psize+1,3) = 0.0D0                 ! v'(r)
    nodes_in(2:psize+1,4) = 0.0D0                 ! v''(r)   
  
    nodes_in(psize+2,1) = end(1)
    nodes_in(psize+2,2) = end(2)
    nodes_in(psize+2,3) = end(3)
    nodes_in(psize+2,4) = 0.0D0


    ! IF EXACTLY A NODE
    complete = .FALSE.
    DO n = 1, SIZE(nodes_in, 1)     
    	IF(r .EQ. nodes_in(n,1))THEN      
  	  	complete = .TRUE.
  	  	y = nodes_in(n,2)
    	END IF
    END DO
 

    ! IF NOT A NODE FIND 
    IF(complete .EQV. .FALSE.)THEN

    	! INTERPOLATE DERIVATIVE INPUT NODES (EXCEPT ZBL START AND FIXED NODE AT END)
    	DO n = 2, SIZE(nodes_in, 1) - 1
    		j = MAX(na,MIN(n,SIZE(nodes_in, 1)-(interp_size-1)))  
    		CALL interpndydxn(nodes_in(n, 1), nodes_in(j:j+(interp_size-1), 1), &
  				 nodes_in(j:j+(interp_size-1), 2), 1, nodes_in(n, 3))
    	END DO 

    	! Option to increase size of nodes, by 4 point interpolation, increase number of nodes for spline   
    	na = 1
    	nb = nsize * (SIZE(nodes_in, 1) - 1) + 1
    	nodes(na:nb,:) = 0.0D0       

    	! INTERPOLATE IF HIGHER LEVEL REQUIRED
    	IF(nsize .EQ. 1)THEN
    		nb = SIZE(nodes_in, 1)
    		nodes(1:nb,:) = nodes_in(:,:)
    	ELSE
    		n = 0
    		DO i = 2, SIZE(nodes_in, 1) - 1
    			DO k = 1, nsize
    				n = n + 1
    				IF(k .EQ. 1)THEN
  	  				nodes(n,:) = nodes_in(i,:)
  	  			ELSE          
  	  				nodes(n,1) = nodes_in(i,1)  + ((k - 1.0D0)/ (1.0D0 * nsize)) * (nodes_in(i+1,1) - nodes_in(i,1))
  	  				j = MAX(1,MIN(i,SIZE(nodes_in, 1)-(interp_size-1)))    
  	  				CALL interpn(nodes(n, 1), nodes_in(j:j+(interp_size-1), 1), &
  	  						 nodes_in(j:j+(interp_size-1), 2), nodes(n, 2))  
  		  			CALL interpn(nodes(n, 1), nodes_in(j:j+(interp_size-1), 1), &
    							 nodes_in(j:j+(interp_size-1), 3), nodes(n, 3))  
    				END IF
    			END DO  
    		END DO
    		n = n + 1
    		nodes(n,:) = nodes_in(SIZE(nodes_in, 1),:)
    		nb = n
    	END IF

    	! FIND n, n+1 TO SPLINE BETWEEN
    	IF(r .LT. nodes(1,1))THEN
    		n = 1    
    	ELSE IF(r .GT. nodes(nb, 1))THEN
    		n = nb  - 1  
    	ELSE
    		DO n = 1, nb - 1
    			IF(r .GT. nodes(n, 1) .AND. r .LT. nodes(n+1, 1))THEN
    				EXIT
  	  		END IF
    		END DO
    	END IF  

    	! X VALUES IN p_fixed
    	! Y VALUES IN p
    	! DERIV VALUES IN dydr
  	
    	xmat(1,1) = 1.0
    	xmat(1,2) = nodes(n, 1)
    	xmat(1,3) = nodes(n, 1)**2
    	xmat(1,4) = nodes(n, 1)**3
    	xmat(2,1) = 1.0
    	xmat(2,2) = nodes(n+1, 1)
    	xmat(2,3) = nodes(n+1, 1)**2
    	xmat(2,4) = nodes(n+1, 1)**3
    	xmat(3,1) = 0.0
    	xmat(3,2) = 1.0D0
    	xmat(3,3) = 2.0D0 * nodes(n, 1)
    	xmat(3,4) = 3.0D0 * nodes(n, 1)**2
    	xmat(4,1) = 0.0
    	xmat(4,2) = 1.0D0
    	xmat(4,3) = 2.0D0 * nodes(n+1, 1)
    	xmat(4,4) = 3.0D0 * nodes(n+1, 1)**2
  		
    	ymat(1) = nodes(n, 2)
    	ymat(2) = nodes(n+1, 2)
    	ymat(3) = nodes(n, 3)
    	ymat(4) = nodes(n+1, 3)
  		
    	CALL sls_solve(xmat, ymat, c)
    	y = c(1) + c(2) * r + c(3) * r**2 + c(4) * r**3

    END IF
  END IF
  
END IF

END SUBROUTINE cubic_knot_spline_fixed_end_pair


! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_fixed_end_pair_v(r, p, p_fixed, y)
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
  CALL cubic_knot_spline_fixed_end_pair(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_fixed_end_pair_v















