SUBROUTINE energy()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################

!$OMP PARALLEL
!$OMP DO
DO cn = 1, cc
  CALL energy_inner(cn)
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
ALLOCATE(density(1:cc, 1:group_max))
config_energy(cn, :) = 0.0D0
density = 0.0D0
DO n = 1, nc
  ! PAIR ENERGY
  CALL pot_search_a(1, nl_l(n, 2), nl_l(n, 4), nl_r(n, 1), fx)
  config_energy(cn, 1) = config_energy(cn, 1) + fx
  
  
  ! A DENSITY AT B
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 2), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 2), f)
    CALL pot_search_a(2, nl_l(n, 2), fn, nl_r(n, 1), fx)  
    density(nl_l(n, 3), fn) = density(nl_l(n, 3), fn) + fx    
    f = f + 1
  END DO 
  
  ! B DENSITY AT A
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 4), f)
    CALL pot_search_a(2, nl_l(n, 4), fn, nl_r(n, 1), fx)
    density(nl_l(n, 1), fn) = density(nl_l(n, 1), fn) + fx
    f = f + 1
  END DO 
END DO
! EMBED ENERGY
DO n = 1, cc
  f = 1
  DO WHILE(fgroups_embe(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 2), f)
    CALL pot_search_a(3, c_l(n), fn, density(n, fn), fx)
    config_energy(cn, 2) = config_energy(cn, 2) + fx
    f = f + 1
  END DO
END DO
DEALLOCATE(density)
config_energy(cn, 3) = config_energy(cn, 1) + config_energy(cn, 2)
config_energy(cn, 4:6) = config_energy(cn, 1:3) / cc
END SUBROUTINE energy_calc
