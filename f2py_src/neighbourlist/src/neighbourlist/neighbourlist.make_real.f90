SUBROUTINE make_real()
!###########################################################
INTEGER(kind=StandardInteger) :: cn, i
!###########################################################

DO cn = 1, cc
  DO i = coords_key(cn,1), coords_key(cn,2)
    coords_real(i, 1:3) = a0(cn) * matmul(uv(cn,1:3,1:3), coords_crystal(i, 1:3)) 
  END DO
END DO


END SUBROUTINE make_real