!##########################################################
! Refresh crystal coords


SUBROUTINE refresh_coords()
!###########################################################
! MAKE CARTESIAN 
INTEGER(KIND=StandardInteger) :: n, i
!###########################################################
DO n = 1, atom_count
  !DO i = 1,3
  !  CALL Modulus(coords(n, i), 1.0D0, coords(n, i))
  !END DO
  coords(n, 1:3) = matmul(from_cartesian, coords(n, 4:6))
END DO  
END SUBROUTINE refresh_coords