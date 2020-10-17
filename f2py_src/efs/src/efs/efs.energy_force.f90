SUBROUTINE energy_force()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################

CALL start_t()
!$OMP PARALLEL
!$OMP DO
DO cn = 1, cc
  CALL energy_force_inner(cn)
END DO
!$OMP END DO
!$OMP END PARALLEL
CALL end_t(ef_timer)
ef_timer_sum = ef_timer_sum + ef_timer
END SUBROUTINE energy_force



SUBROUTINE energy_force_inner(cn)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################
CALL energy_force_calc(cn, &
                 nlist_l(key(cn, 5):key(cn, 6),:), &
                 nlist_r(key(cn, 5):key(cn, 6),:), &
                 labels(key(cn, 1):key(cn, 2)), &  
                 coords(key(cn, 1):key(cn, 2), :), &  
                 config_forces(key(cn, 1):key(cn, 2), :) & 
                )
END SUBROUTINE energy_force_inner


SUBROUTINE energy_force_calc(cn, nl_l, nl_r, c_l, c_r, efs_forces)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
INTEGER(KIND=StandardInteger) :: nl_l(:,:)
REAL(kind=DoubleReal) :: nl_r(:,:)
INTEGER(KIND=StandardInteger) :: c_l(:)
REAL(kind=DoubleReal) :: c_r(:,:)
REAL(kind=DoubleReal) :: efs_forces(:,:)
!###########################################################
INTEGER(KIND=StandardInteger) :: n, f, fn, dn
INTEGER(KIND=StandardInteger) :: cc, nc
REAL(kind=DoubleReal) :: fx(1:2)
REAL(kind=DoubleReal) :: fv(1:3)
REAL(kind=DoubleReal), ALLOCATABLE :: density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ab(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ba(:,:)
REAL(kind=DoubleReal) :: epA, epB, dpAB, dpBA
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
ALLOCATE(density_grad_ab(1:nc, 1:group_max))
ALLOCATE(density_grad_ba(1:nc, 1:group_max))
config_energy(cn, :) = 0.0D0
efs_forces(:,:) = 0.0D0
density = 0.0D0
density_grad_ab = 0.0D0
density_grad_ba = 0.0D0
DO n = 1, nc
  ! PAIR ENERGY
  CALL pot_search_b(1, nl_l(n, 2), nl_l(n, 4), nl_r(n, 1), fx)
  config_energy(cn, 1) = config_energy(cn, 1) + fx(1)
  
  ! PAIR FORCE
  fv(:) = fx(2) * nl_r(n, 2:4)
  efs_forces(nl_l(n, 1), :) = efs_forces(nl_l(n, 1), :) - fv(:)
  efs_forces(nl_l(n, 3), :) = efs_forces(nl_l(n, 3), :) + fv(:)
  
  
  ! A DENSITY AT B
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 2), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 2), f)
    CALL pot_search_b(2, nl_l(n, 2), fn, nl_r(n, 1), fx)  
    density(nl_l(n, 3), fn) = density(nl_l(n, 3), fn) + fx(1) 
    density_grad_ab(n, fn) = density_grad_ab(n, fn) + fx(2) 
    f = f + 1
  END DO 
  
  ! B DENSITY AT A
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 4), f)
    CALL pot_search_b(2, nl_l(n, 4), fn, nl_r(n, 1), fx)
    density(nl_l(n, 1), fn) = density(nl_l(n, 1), fn) + fx(1)
    density_grad_ba(n, fn) = density_grad_ba(n, fn) + fx(2) 
    f = f + 1
  END DO 
END DO

! EMBED ENERGY
DO n = 1, cc
  f = 1
  DO WHILE(fgroups_embe(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_embe(nl_l(n, 2), f)
    CALL pot_search_b(3, c_l(n), fn, density(n, fn), fx)
    config_energy(cn, 2) = config_energy(cn, 2) + fx(1)
    f = f + 1    
  END DO
END DO



DO n = 1, nc
  f = 1
  DO WHILE(fgroups_embe(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_embe(nl_l(n, 2), f)
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_b(3, nl_l(n, 2), fn, density(nl_l(n, 1), fn), fx)
    epA = fx(2)
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_b(3, nl_l(n, 4), fn, density(nl_l(n, 3), fn), fx)
    epB = fx(2)
    
    ! @Pij(r)/@r   Gradient of density function at A due to atom B
    dpAB = density_grad_ab(n,fn)
      
    ! @Pji(r)/@r   Gradient of density function at B due to atom A
    dpBA = density_grad_ba(n,fn)

  
    fv(:) = (epA * dpBA + epB * dpAB) * nl_r(n, 2:4)
    efs_forces(nl_l(n, 1), :) = efs_forces(nl_l(n, 1), :) - fv(:)
    efs_forces(nl_l(n, 3), :) = efs_forces(nl_l(n, 3), :) + fv(:)
  
    f = f + 1    
  END DO
END DO


DEALLOCATE(density)
DEALLOCATE(density_grad_ab)
DEALLOCATE(density_grad_ba)
config_energy(cn, 3) = config_energy(cn, 1) + config_energy(cn, 2)
config_energy(cn, 4:6) = config_energy(cn, 1:3) / cc
END SUBROUTINE energy_force_calc
