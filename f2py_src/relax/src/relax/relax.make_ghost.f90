!##########################################################
! Makes ghost cell


SUBROUTINE make_ghost()
!###########################################################
INTEGER(kind=StandardInteger) :: ccz
REAL(kind=DoubleReal) :: to_cart(1:3, 1:3)

INTEGER(KIND=StandardInteger) :: cn, an, cx, cy, cz, n, m, i

REAL(kind=DoubleReal) :: shift(1:3)
REAL(kind=DoubleReal) :: coords_crystal(1:3)
REAL(kind=DoubleReal) :: coords_cartesian(1:3)
REAL(kind=DoubleReal) :: half(1:3)
REAL(kind=DoubleReal) :: c(1:3)
REAL(kind=DoubleReal) :: ct(1:3)
REAL(kind=DoubleReal) :: td(1:3)
REAL(kind=DoubleReal) :: rd(1:3) 
REAL(kind=DoubleReal) :: r(1:3) 
REAL(kind=DoubleReal) :: vol_cell(1:3,1:3)
!###########################################################

half = 0.5D0
c = alat * MATMUL(half, uv(1:3 , 1:3))
td = gthreshold * (c + rcut)
    
an = 1
DO cx = -1, 1
  DO cy = -1 ,1
    DO cz = -1, 1
      IF(cx .EQ. 0 .AND. cy .EQ. 0 .AND. cz .EQ. 0) THEN
        DO n = 1, atom_count
          ghostn(an) = n
          ghostids(an) = ids(n)
          ghostlabels(an) = labels(n)
          ghostcoords(an, 1:3) = coords(n, 4:6)
          ghostshift(an, 1:3) = 0.0D0
          ghosthalo(an) = .FALSE. 
          an = an + 1
        END DO
      ELSE     
        ! Shift
        shift(1) = 1.0D0 * cx
        shift(2) = 1.0D0 * cy
        shift(3) = 1.0D0 * cz
        
        DO n = 1, atom_count 
          ct(1:3) = matmul(alat * uv(1:3, 1:3), coords(n, 1:3) + shift(1:3))
          rd(1:3) = abs(ct(1:3) - c(1:3))          
          IF(rd(1) .LE. td(1) .AND. rd(2) .LE. td(2) .AND. rd(3) .LE. td(3))THEN
            ghostn(an) = n
            ghostids(an) = ids(n)
            ghostlabels(an) = labels(n)
            ghostcoords(an, 1:3) = ct(1:3)
            ghostshift(an, 1:3) = shift(1:3)
            ghosthalo(an) = .TRUE.  
            an = an + 1
          END IF
        END DO 
      END IF   
    END DO
  END DO
END DO
ghost_count = an - 1

END SUBROUTINE make_ghost

