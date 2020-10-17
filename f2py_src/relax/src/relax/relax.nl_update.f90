SUBROUTINE nl_update()
!###########################################################
REAL(kind=DoubleReal) :: r_cut_sq, r_sq
REAL(kind=DoubleReal) :: r(1:3), rd(1:3), c(1:3), r_mag
INTEGER(KIND=StandardInteger) :: an, gn, n, nc
!###########################################################

DO nc = 1, nl_count
  r(1:3) = ghostcoords(nlist_l(nc, 6), 1:3) - coords(nlist_l(nc, 5), 4:6)
  r_mag = sqrt(SUM(r(1:3) * r(1:3)))     
  nlist_r(nc, 1) = r_mag
  nlist_r(nc, 2:4) = r(1:3)/r_mag
END DO


END SUBROUTINE nl_update
