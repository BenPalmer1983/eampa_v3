SUBROUTINE get_energy(cn, energy_out)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: cn
REAL(kind=DoubleReal), INTENT(OUT) :: energy_out
!###########################################################
energy_out = config_energy(cn, 3)
END SUBROUTINE get_energy



SUBROUTINE get_energies(energy_out)
!###########################################################
REAL(kind=DoubleReal), INTENT(OUT) :: energy_out(1:1000)
!###########################################################
energy_out(:) = config_energy(:,3)
END SUBROUTINE get_energies


FUNCTION get_force(cn) RESULT (force_out)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
REAL(kind=DoubleReal), ALLOCATABLE :: force_out(:,:)
!###########################################################
IF(ALLOCATED(force_out))THEN
  DEALLOCATE(force_out)
END IF
print *, "A"
ALLOCATE(force_out(1:key(cn, 2)-key(cn, 1)+1, 1:3))
force_out(:,:) = config_forces(key(cn, 1):key(cn, 2), :)
print *, force_out
END FUNCTION get_force