SUBROUTINE energy()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################
max_density = 0.0D0
!$OMP PARALLEL
!$OMP DO
DO cn = 1, cc
  CALL energy_inner(cn)
  !print *, cn, config_energy(cn, 3)
  !print *, config_energy(cn, 1), config_energy(cn, 2), config_energy(cn, 3), config_energy(cn, 6)
END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE energy



SUBROUTINE energy_inner(cn)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################
CALL energy_calc(cn, &
                 nlist_l(key(cn, 5):key(cn, 6),:), &
                 nlist_r(key(cn, 5):key(cn, 6),:), &
                 labels(key(cn, 1):key(cn, 2)), &  
                 coords(key(cn, 1):key(cn, 2), :) & 
                )
END SUBROUTINE energy_inner


SUBROUTINE energy_calc(cn, nl_l, nl_r, c_l, c_r)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
INTEGER(KIND=StandardInteger) :: nl_l(:,:)
REAL(kind=DoubleReal) :: nl_r(:,:)
INTEGER(KIND=StandardInteger) :: c_l(:)
REAL(kind=DoubleReal) :: c_r(:,:)
!###########################################################
INTEGER(KIND=StandardInteger) :: n, f, fn, dn
INTEGER(KIND=StandardInteger) :: cc, nc
REAL(kind=DoubleReal) :: fx
REAL(kind=DoubleReal), ALLOCATABLE :: density(:,:)
!###########################################################
! nlist_l(nc, 1) = ids(an)
! nlist_l(nc, 2) = labels(an)
! nlist_l(nc, 3) = ghostids(gn)
! nlist_l(nc, 4) = ghostlabels(gn)          
! nlist_r(nc, 1) = r_mag
! nlist_r(nc, 2:4) = r(1:3)/r_mag
cc = SIZE(c_l, 1)
nc = SIZE(nl_l, 1)
ALLOCATE(density(1:cc, 1:fgroup_max))
config_energy(cn, :) = 0.0D0
density = 0.0D0
DO n = 1, nc
  ! PAIR ENERGY
  CALL search_f(1, nl_l(n, 2), nl_l(n, 4), 0, nl_r(n, 1), fx)
  config_energy(cn, 1) = config_energy(cn, 1) + fx
  ! LOOP THROUGH DENSITY GROUPS
  DO fn = 1, fgroup_max
    ! DENSITY OF ATOM B AT A (B with any atom)
    CALL search_f(2, nl_l(n, 2), 0, fn, nl_r(n, 1), fx)
    density(nl_l(n, 1), fn) = density(nl_l(n, 1), fn) + fx 
    ! DENSITY OF ATOM B AT A (B with only atom A)
    CALL search_f(2, nl_l(n, 2), nl_l(n, 4), fn, nl_r(n, 1), fx)
    density(nl_l(n, 1), fn) = density(nl_l(n, 1), fn) + fx 
    ! DENSITY OF ATOM A AT B (A with any atom)
    CALL search_f(2, nl_l(n, 4), 0, fn, nl_r(n, 1), fx)
    density(nl_l(n, 3), fn) = density(nl_l(n, 3), fn) + fx 
    ! DENSITY OF ATOM A AT B (A with only atom B)
    CALL search_f(2, nl_l(n, 4), nl_l(n, 2), fn, nl_r(n, 1), fx)
    density(nl_l(n, 3), fn) = density(nl_l(n, 3), fn) + fx    
  END DO
END DO
! EMBED ENERGY
DO n = 1, cc
  ! LOOP THROUGH DENSITY GROUPS
  DO fn = 1, fgroup_max    
    CALL search_f(3, nl_l(n, 4), 0, fn, density(n, fn), fx)
    config_energy(cn, 2) = config_energy(cn, 2) + fx
  END DO
END DO
max_density = MAX(max_density, MAXVAL(density))
DEALLOCATE(density)
config_energy(cn, 3) = config_energy(cn, 1) + config_energy(cn, 2)
config_energy(cn, 4:6) = config_energy(cn, 1:3) / cc

END SUBROUTINE energy_calc
