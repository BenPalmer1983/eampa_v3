SUBROUTINE energy_force_stress()
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################

CALL start_t()
!$OMP PARALLEL
!$OMP DO
DO cn = 1, cc
  CALL energy_inner(cn)
  CALL energy_force_stress_inner(cn)
END DO
!$OMP END DO
!$OMP END PARALLEL
CALL end_t(efs_timer)
efs_timer_sum = efs_timer_sum + efs_timer
END SUBROUTINE energy_force_stress



SUBROUTINE energy_force_stress_inner(cn)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
!###########################################################
CALL energy_force_stress_calc(cn, &
                 nlist_l(key(cn, 5):key(cn, 6),:), &
                 nlist_r(key(cn, 5):key(cn, 6),:), &
                 labels(key(cn, 1):key(cn, 2)), &  
                 coords(key(cn, 1):key(cn, 2),:) &
                )
END SUBROUTINE energy_force_stress_inner


SUBROUTINE energy_force_stress_calc(cn, nl_l, nl_r, c_l, c_r)
!###########################################################
INTEGER(KIND=StandardInteger) :: cn
INTEGER(KIND=StandardInteger) :: nl_l(:,:)
REAL(kind=DoubleReal) :: nl_r(:,:)
INTEGER(KIND=StandardInteger) :: c_l(:)
REAL(kind=DoubleReal) :: c_r(:,:)
!###########################################################
INTEGER(KIND=StandardInteger) :: i, j, n, f, fn, dn
INTEGER(KIND=StandardInteger) :: cc, nc
REAL(kind=DoubleReal) :: fx
REAL(kind=DoubleReal) :: fx_arr(1:2)
REAL(kind=DoubleReal) :: fv(1:3)
REAL(kind=DoubleReal), ALLOCATABLE :: density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ab(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ba(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: nl_force(:,:)
!REAL(kind=DoubleReal) :: s(1:3, 1:3)
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
ALLOCATE(density(1:cc, 1:fgroup_max))
ALLOCATE(density_grad_ab(1:nc, 1:fgroup_max))
ALLOCATE(density_grad_ba(1:nc, 1:fgroup_max))
ALLOCATE(nl_force(1:nc, 1:3))
config_energy(cn, :) = 0.0D0
config_forces(cn, 1:key(cn, 20), 1:3) = 0.0D0
density = 0.0D0
density_grad_ab = 0.0D0
density_grad_ba = 0.0D0

DO n = 1, nc
  ! PAIR ENERGY
  CALL pot_search_pair_arr(nl_l(n, 2), nl_l(n, 4), nl_r(n, 1), fx_arr)
  config_energy(cn, 1) = config_energy(cn, 1) + fx_arr(1)
  
  ! PAIR FORCE
  fv(:) = fx_arr(2) * nl_r(n, 2:4)
  config_forces(cn, nl_l(n, 1), :) = config_forces(cn, nl_l(n, 1), :) - fv(:)
  config_forces(cn, nl_l(n, 3), :) = config_forces(cn, nl_l(n, 3), :) + fv(:)
  nl_force(n,:) = - fv(:)
  
  ! A DENSITY AT B
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 2), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 2), f)
    CALL pot_search_dens_arr(nl_l(n, 2), fn, nl_r(n, 1), fx_arr)  
    density(nl_l(n, 3), fn) = density(nl_l(n, 3), fn) + fx_arr(1) 
    density_grad_ab(n, fn) = density_grad_ab(n, fn) + fx_arr(2) 
    f = f + 1
  END DO 
  
  ! B DENSITY AT A
  f = 1
  DO WHILE(fgroups_dens(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_dens(nl_l(n, 4), f)
    CALL pot_search_dens_arr(nl_l(n, 4), fn, nl_r(n, 1), fx_arr)
    density(nl_l(n, 1), fn) = density(nl_l(n, 1), fn) + fx_arr(1)
    density_grad_ba(n, fn) = density_grad_ba(n, fn) + fx_arr(2) 
    f = f + 1
  END DO 
END DO

f = 1
DO WHILE(fgroups_embe(nl_l(1, 4), f) .NE. 0)
  fn = fgroups_embe(nl_l(1, 2), f)
  f = f + 1    
END DO


! EMBED ENERGY
DO n = 1, cc
  f = 1
  DO WHILE(fgroups_embe(nl_l(n, 4), f) .NE. 0)
    fn = fgroups_embe(nl_l(n, 2), f)
    CALL pot_search_embe(c_l(n), fn, density(n, fn), fx)
    config_energy(cn, 2) = config_energy(cn, 2) + fx
    f = f + 1    
  END DO
END DO


! EMBED FORCE
DO n = 1, nc
  f = 1
  DO WHILE(fgroups_embe(nl_l(n, 4), f) .NE. 0)

    fn = fgroups_embe(nl_l(n, 2), f)
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_embe_arr(nl_l(n, 2), fn, density(nl_l(n, 1), fn), fx_arr)
    epA = fx_arr(2)
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_embe_arr(nl_l(n, 4), fn, density(nl_l(n, 3), fn), fx_arr)
    epB = fx_arr(2)
    
    ! @Pij(r)/@r   Gradient of density function at A due to atom B
    dpAB = density_grad_ab(n,fn)
      
    ! @Pji(r)/@r   Gradient of density function at B due to atom A
    dpBA = density_grad_ba(n,fn)

  
    fv(:) = (epA * dpBA + epB * dpAB) * nl_r(n, 2:4)
    config_forces(cn, nl_l(n, 1), :) = config_forces(cn, nl_l(n, 1), :) - fv(:)
    config_forces(cn, nl_l(n, 3), :) = config_forces(cn, nl_l(n, 3), :) + fv(:)
  
    nl_force(n,:) = nl_force(n,:) - fv(:)
    f = f + 1    
  END DO
END DO


! STRESS
config_stresses(cn, 1:3, 1:3) = 0.0D0
DO n = 1, nc
  IF(nlisthalo(n))THEN  
    DO i = 1,3
      DO j = 1,3
        config_stresses(cn, i,j) = config_stresses(cn, i,j) + (nl_r(n, 1+i) * nl_force(n,j))
      END DO
    END DO
  END IF
END DO
config_stresses(cn, 1:3, 1:3) = config_stresses(cn, 1:3, 1:3) / (2.0D0 * volume(cn))

max_density = MAX(max_density, MAXVAL(density))


DEALLOCATE(density)
DEALLOCATE(density_grad_ab)
DEALLOCATE(density_grad_ba)
DEALLOCATE(nl_force)
config_energy(cn, 3) = config_energy(cn, 1) + config_energy(cn, 2)
config_energy(cn, 4:6) = config_energy(cn, 1:3) / cc

!print *, cn
!print *, config_energy(cn, 1)
!print *, config_energy(cn, 2)
!print *, config_energy(cn, 3)

END SUBROUTINE energy_force_stress_calc


! OLD STRESS
!DO n = 1, nc
!  DO i = 1,3
!    DO j = 1,3
!      s(i,j) = s(i,j) + (nl_r(n, 1+i) * nl_force(n,j))
!      IF(.NOT. nlisthalo(n))THEN        
!        s(i,j) = s(i,j) - (nl_r(n, 1+i) * nl_force(n,j))
!      END IF
!    END DO
!  END DO
!END DO