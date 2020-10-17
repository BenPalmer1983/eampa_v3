!##########################################################
! Makes ghost cell


SUBROUTINE ghost_update()
!###########################################################
INTEGER(kind=StandardInteger) :: gn, n
!###########################################################

DO gn = 1, ghost_count
  n = ghostn(gn)
  ghostcoords(gn, 1:3) = matmul(alat * uv(1:3, 1:3), coords(n, 1:3) + ghostshift(gn, 1:3))
END DO

END SUBROUTINE ghost_update

