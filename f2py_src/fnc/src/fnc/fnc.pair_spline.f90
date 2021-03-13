! PAIR SPLINE
! Ackland-Mendelev ZBL directly into a spline
! 
! Pfixed    q1, q2, zbl_cutoff  node1 node 2...node n end_node

! SCALAR SUBROUTINE
SUBROUTINE pair_spline(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: Q1, Q2, zbl_r, end_node_r
REAL(kind=DoubleReal) :: rs, e, x, ep, epp
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb, H, zbl_dydr
REAL(kind=DoubleReal) :: start_x, start_y, start_yp
REAL(kind=DoubleReal) :: end_x, end_y, end_yp
INTEGER(kind=StandardInteger) :: i, j, n, n_p, n_pf, interp_size
REAL(kind=DoubleReal) :: nodes(1:SIZE(p_fixed,1) - 2, 1:3)
LOGICAL :: complete
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
!############################################################
Q1 = p_fixed(SIZE(p_fixed,1)-3)
Q2 = p_fixed(SIZE(p_fixed,1)-2)
zbl_r = p_fixed(SIZE(p_fixed,1)-1)
end_node_r = p_fixed(SIZE(p_fixed,1))
interp_size = 2*(SIZE(p_fixed,1))

IF(SIZE(p_fixed,1) .NE. SIZE(p,1) + 4)THEN
  y = 0.0D0
ELSE IF(r .EQ. 0.0D0)THEN
  y = 1.0D4
ELSE IF(r .GT. 0.0D0 .AND. r .LE. zbl_r)THEN
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (r / rs)
  e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
  y = ((Q1 * Q2) / r) * e
ELSE IF(r .GT. zbl_r .AND. r .LT. end_node_r)THEN  
  ! INTERPOLATE INNER NODES DY/DR
  nodes = 0.0D0
  
  ! ZBL Node
  start_x = zbl_r
  rs = 0.4683766D0 / (Q1**(2.0/3.0D0) + Q2**(2.0/3.0D0))
  x = (zbl_r / rs)
  e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
  start_y = ((Q1 * Q2) / zbl_r) * e
  rs = 0.4683766 / (Q1**(2.0D0/3.0D0)+Q2**(2.0D0/3.0D0))   
  ep =      0.1818D0 * exp(zbl_r * (-3.2D0/rs)) 
  ep = ep + 0.5099D0 * exp(zbl_r * (-0.9423D0/rs)) 
  ep = ep + 0.2802 * exp(zbl_r * (-0.4029D0/rs)) 
  ep = ep + 0.02817 * exp(zbl_r * (-0.2016D0/rs))     
  epp =      (-3.2D0/rs) * 0.1818D0 * exp(zbl_r * (-3.2D0/rs)) 
  epp = ep + (-0.9423D0/rs) * 0.5099D0 * exp(zbl_r * (-0.9423D0/rs)) 
  epp = ep + (-0.4029D0/rs) * 0.2802 * exp(zbl_r * (-0.4029D0/rs)) 
  epp = ep + (-0.2016D0/rs) * 0.02817 * exp(zbl_r * (-0.2016D0/rs))
  start_yp = ((Q1 * Q2) / zbl_r) * epp - ((Q1 * Q2) / (zbl_r**2)) * ep
  
  ! END NODE
  end_x = end_node_r
  end_y = 0.0D0
  end_yp = 0.0D0
  
  ! Start Node
  nodes(1, 1) = start_x
  nodes(1, 2) = start_y
  nodes(1, 3) = start_yp

  ! Inner nodes
  DO n = 1, SIZE(nodes, 1)-2
    nodes(n+1, 1) = p_fixed(n)
  END DO  
  DO n = 1, SIZE(nodes, 1)-2
    nodes(n+1, 2) = p(n)
  END DO

  ! End node
  nodes(SIZE(nodes, 1), 1) = end_x
  nodes(SIZE(nodes, 1), 2) = end_y
  nodes(SIZE(nodes, 1), 3) = end_yp
  
  DO n = 2, SIZE(nodes, 1)-1
    j = n-2
    IF(j .LT. 1)THEN
      j = 1
    ELSE IF(j+3 .GT. SIZE(nodes, 1))THEN
      j = SIZE(nodes, 1) - 3
    END IF
    CALL interp4dydx(nodes(n,1), nodes(j:j+3,1), nodes(j:j+3,2), nodes(n,3))
  END DO
  

  ! IF EXACTLY A NODE
  complete = .FALSE.
  DO n = 1, SIZE(nodes, 1)  
    IF(r .EQ. nodes(n, 1))THEN      
      complete = .TRUE.
      y = nodes(n, 2) 
      EXIT
    END IF
  END DO
 
  ! IF BETWEEN NODES
  IF(complete .EQV. .FALSE.)THEN  
    
    DO n = 1, SIZE(nodes, 1)-1
      IF(r .GT. nodes(n,1) .AND. r .LT. nodes(n+1,1))THEN
        EXIT
      END IF
    END DO
    
    xmat(1,1) = 1.0
    xmat(1,2) = nodes(n,1)
    xmat(1,3) = nodes(n,1)**2
    xmat(1,4) = nodes(n,1)**3
    xmat(2,1) = 1.0
    xmat(2,2) = nodes(n+1,1)
    xmat(2,3) = nodes(n+1,1)**2
    xmat(2,4) = nodes(n+1,1)**3
    xmat(3,1) = 0.0
    xmat(3,2) = 1.0D0
    xmat(3,3) = 2.0D0 * nodes(n,1)
    xmat(3,4) = 3.0D0 * nodes(n,1)**2
    xmat(4,1) = 0.0
    xmat(4,2) = 1.0D0
    xmat(4,3) = 2.0D0 * nodes(n+1,1)
    xmat(4,4) = 3.0D0 * nodes(n+1,1)**2
      
    ymat(1) = nodes(n,2)
    ymat(2) = nodes(n+1,2)
    ymat(3) = nodes(n,3)
    ymat(4) = nodes(n+1,3)
      
    CALL sls_solve(xmat, ymat, c)
    y = c(1) + c(2) * r + c(3) * r**2 + c(4) * r**3

  END IF

ELSE IF(r .EQ. end_node_r)THEN
  y = 0.0D0
ELSE
  y = 0.0D0
END IF


END SUBROUTINE pair_spline



! VECTOR SUBROUTINE
SUBROUTINE pair_spline_v(r, p, p_fixed, y)
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
  CALL pair_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE pair_spline_v




SUBROUTINE interp4(xi, x, y, yi)
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
END SUBROUTINE interp4



SUBROUTINE interp4dydx(xi, x, y, ypi)
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
END SUBROUTINE interp4dydx


