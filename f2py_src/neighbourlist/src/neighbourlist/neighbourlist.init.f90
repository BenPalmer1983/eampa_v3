SUBROUTINE init()
!###########################################################

IF(ALLOCATED(atom_id))THEN
  DEALLOCATE(atom_id)
END IF
IF(ALLOCATED(atom_type))THEN
  DEALLOCATE(atom_type)
END IF
IF(ALLOCATED(coords_crystal))THEN
  DEALLOCATE(coords_crystal)
END IF
IF(ALLOCATED(coords_real))THEN
  DEALLOCATE(coords_real)
END IF

IF(ALLOCATED(atom_halo_id))THEN
  DEALLOCATE(atom_halo_id)
END IF
IF(ALLOCATED(atom_halo_type))THEN
  DEALLOCATE(atom_halo_type)
END IF
IF(ALLOCATED(coords_halo))THEN
  DEALLOCATE(coords_halo)
END IF
IF(ALLOCATED(in_halo))THEN
  DEALLOCATE(in_halo)
END IF

IF(ALLOCATED(nl_atoms))THEN
  DEALLOCATE(nl_atoms)
END IF



ALLOCATE(atom_id(1:max_coords))
ALLOCATE(atom_type(1:max_coords))
ALLOCATE(coords_crystal(1:max_coords,1:3))
ALLOCATE(coords_real(1:max_coords,1:3))


ALLOCATE(atom_halo_id(1:max_coords_halo))
ALLOCATE(atom_halo_type(1:max_coords_halo))
ALLOCATE(coords_halo(1:max_coords_halo,1:3))
ALLOCATE(in_halo(1:max_coords_halo))


ALLOCATE(nl_atoms(1:max_nl,1:4))
ALLOCATE(nl_r(1:max_nl))
ALLOCATE(nl_rvec(1:max_nl,1:3))
ALLOCATE(nl_in_halo(1:max_nl))






END SUBROUTINE init