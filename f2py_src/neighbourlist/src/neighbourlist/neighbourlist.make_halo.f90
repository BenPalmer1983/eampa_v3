SUBROUTINE make_halo()
!###########################################################
INTEGER(kind=StandardInteger) :: cn, i, nx, ny, nz, gn
REAL(kind=DoubleReal) :: shift(1:3)
!###########################################################


gn = 1
DO cn = 1, cc
  halo_key(cn,1) = gn
  DO nx = -1, 1
    DO ny = -1,1
      DO nz = -1,1
        shift(1) = 1.0D0 * nx
	      shift(2) = 1.0D0 * ny
		    shift(3) = 1.0D0 * nz
        DO i = coords_key(cn,1), coords_key(cn,2)
	        IF(nx .EQ. 0 .AND. ny .EQ. 0 .AND. nz .EQ. 0)THEN
	          in_halo(gn) = .FALSE.
	        ELSE
	          in_halo(gn) = .TRUE.
	        END IF
          atom_halo_id(gn) = atom_id(i)
	        atom_halo_type(gn) = atom_type(i)
          coords_halo(gn, 1:3) = a0(cn) * matmul(uv(cn,1:3,1:3), coords_crystal(i, 1:3) + shift(1:3))
          gn = gn + 1
        END DO
      END DO  
    END DO
  END DO  
	halo_key(cn,2) = gn - 1
  halo_key(cn,3) = gn - halo_key(cn,1)
END DO





END SUBROUTINE make_halo