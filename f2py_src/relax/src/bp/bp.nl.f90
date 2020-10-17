SUBROUTINE make_nl(cn)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: cn
!###########################################################
REAL(kind=DoubleReal) :: r_cut_sq, r_sq
REAL(kind=DoubleReal) :: r(1:3), rd(1:3), c(1:3), r_mag
INTEGER(KIND=StandardInteger) :: an, gn, n, nc
!###########################################################

! rcut 
r_cut_sq = rcut(cn) * rcut(cn)

! STORE KEYS
IF(cn .EQ. 1)THEN
  key(cn, 5) = 1
ELSE  
  key(cn, 5) = key(cn-1, 6) + 1
END IF
  
nc = key(cn, 5)  
  
DO an = key(cn, 1), key(cn, 2)
  DO gn = key(cn, 3), key(cn, 4)
    IF (ids(an) .LT. ghostids(gn)) THEN
      r(1:3) = ghostcoords(gn, 1:3) - coords(an, 4:6)
      r_sq = SUM(r(1:3) * r(1:3))
      IF (r_sq .LE. r_cut_sq) THEN          
        r_mag = sqrt(r_sq)          
        nlist_l(nc, 1) = ids(an)
        nlist_l(nc, 2) = labels(an)
        nlist_l(nc, 3) = ghostids(gn)
        nlist_l(nc, 4) = ghostlabels(gn)          
        nlist_r(nc, 1) = r_mag
        nlist_r(nc, 2:4) = r(1:3)/r_mag
        nc = nc + 1
      END IF
    END IF
  END DO
END DO  
key(cn, 6) = nc - 1

END SUBROUTINE make_nl
