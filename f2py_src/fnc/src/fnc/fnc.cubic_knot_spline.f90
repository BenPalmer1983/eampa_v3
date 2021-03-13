! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
! Pfixed    q1, q2, zbl_cutoff  node1 node 2...node n end_node



SUBROUTINE cubic_knot_spline(r, p, p_fixed, y)
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
INTEGER(kind=StandardInteger) :: i, j, k, n, n_p, n_pf, interp_size
REAL(kind=DoubleReal) :: nodes_in(1:SIZE(p,1),1:3)
REAL(kind=DoubleReal) :: nodes(1:10000,1:3)
INTEGER(kind=StandardInteger) :: na, nb, nsize
LOGICAL :: complete
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
!############################################################

! P and PF must be the same length
! PF contains x value of node
! P contains y value of node

IF(SIZE(p_fixed,1) .LT. SIZE(p,1))THEN
  y = 0.0D0
ELSE
  ! LOAD X and Y POINTS
  nodes_in(:,1) = p_fixed(1:SIZE(p,1))
  nodes_in(:,2) = p(:)
  nodes_in(:,3) = 0.0D0
 
  ! IF EXACTLY A NODE
  complete = .FALSE.
  DO n = 1, SIZE(nodes_in, 1)  
    IF(r .EQ. nodes_in(n,1))THEN      
      complete = .TRUE.
      y = p(n) 
      EXIT
    END IF
  END DO

  ! IF NOT A NODE FIND 
  IF(complete .EQV. .FALSE.)THEN
  
    ! Option to increase size of nodes, by 4 point interpolation, increase number of nodes for spline
    nsize = 1
    IF(SIZE(p_fixed,1) .GT. SIZE(p, 1))THEN
      nsize = INT(CEILING(p_fixed(SIZE(p, 1)+1)))
    END IF

    
    na = 1
    nb = nsize * (SIZE(nodes_in, 1) - 1) + 1
    nodes(na:nb,:) = 0.0D0

    ! INTERPOLATE IF HIGHER LEVEL REQUIRED
    IF(nsize .EQ. 1)THEN
      nb = SIZE(nodes_in, 1)
      nodes(1:nb,:) = nodes_in(:,:)
    ELSE
      n = 0
      DO i = 1, SIZE(nodes_in, 1) - 1
        DO k = 1, nsize
          n = n + 1
          IF(k .EQ. 1)THEN
            nodes(n,:) = nodes_in(i,:)
          ELSE          
            nodes(n,1) = nodes_in(i,1)  + ((k - 1.0D0)/ (1.0D0 * nsize)) * (nodes_in(i+1,1) - nodes_in(i,1))
            j = MAX(1,MIN(i,SIZE(nodes_in, 1)-3))    
            CALL cubic_knot_spline_interp4(nodes(n, 1), nodes_in(j:j+3, 1), nodes_in(j:j+3, 2), nodes(n, 2))
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
    ! FIND 4 NEAREST NODES TO INTERPOLATE DERIVATIVE
    j = n - 1
    IF(j .LT. 1)THEN
      j = 1
    ELSE IF(j+3 .GT. nb)THEN
      j = nb - 3
    END IF
    DO k=j,j+3
      CALL cubic_knot_spline_interp4dydx(nodes(k, 1), nodes(j:j+3, 1), nodes(j:j+3, 2), nodes(k, 3))
    END DO

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

END SUBROUTINE cubic_knot_spline



! VECTOR SUBROUTINE
SUBROUTINE cubic_knot_spline_v(r, p, p_fixed, y)
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
  CALL cubic_knot_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE cubic_knot_spline_v




SUBROUTINE cubic_knot_spline_interp4(xi, x, y, yi)
! Identity for square matrix
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: y(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: yi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j
REAL(kind=DoubleReal) :: li
!############################################################
yi = 0.0D0
IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
  DO i = 1, 4
    li = 1.0D0
    DO j = 1, 4
      IF(i .NE. j) THEN
        li = li * (xi - x(j)) / (x(i) - x(j))
      END IF
    END DO
    yi = yi + li * y(i)
  END DO
END IF
!############################################################
END SUBROUTINE cubic_knot_spline_interp4



SUBROUTINE cubic_knot_spline_interp4dydx(xi, x, y, ypi)
! Interpolate and return derivative at xi
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: y(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: ypi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, k, n
REAL(kind=DoubleReal) :: fx, gx, psum
!############################################################
n = 4
IF (SIZE(x,1) .EQ. SIZE(y,1) .AND. n .EQ. SIZE(x,1)) THEN
  ypi = 0.0D0
  Do i=1,SIZE(x,1)
    fx = 1.0D0
    gx = 0.0D0
    Do j=1,SIZE(x,1)
      If(i .NE. j) Then
        fx = fx / (x(i) - x(j))
        psum = 1.0D0
        Do k=1,SIZE(x,1)
          If((i .NE. k) .AND. (j .NE. k))Then
            psum = psum * (xi - x(k))
          End If
        End Do
        gx = gx + psum
      End If
    End Do
    ypi = ypi + fx * gx * y(i)
  End Do
END IF
!############################################################
END SUBROUTINE cubic_knot_spline_interp4dydx


