SUBROUTINE force_calc()
!###########################################################
INTEGER(KIND=StandardInteger) :: n, f, fn, dn
REAL(kind=DoubleReal) :: fx
REAL(kind=DoubleReal) :: fx_d(1:2)
REAL(kind=DoubleReal) :: fv(1:3)
REAL(kind=DoubleReal), ALLOCATABLE :: density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ab(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: density_grad_ba(:,:)
REAL(kind=DoubleReal) :: epA, epB, dpAB, dpBA, e
!###########################################################
! nlist_l(nc, 1) = ids(an)
! nlist_l(nc, 2) = labels(an)
! nlist_l(nc, 3) = ghostids(gn)
! nlist_l(nc, 4) = ghostlabels(gn)          
! nlist_r(nc, 1) = r_mag
! nlist_r(nc, 2:4) = r(1:3)/r_mag


ALLOCATE(density(1:atom_count, 1:group_max))
ALLOCATE(density_grad_ab(1:nl_count, 1:group_max))
ALLOCATE(density_grad_ba(1:nl_count, 1:group_max))
forces(1:atom_count,:) = 0.0D0
density = 0.0D0
density_grad_ab = 0.0D0
density_grad_ba = 0.0D0
e = 0.0D0

DO n = 1, nl_count
  ! PAIR FORCE
  CALL pot_search_a(1, nlist_l(n, 2), nlist_l(n, 4), nlist_r(n, 1), fx)
  e = e + fx
  CALL pot_search_c(1, nlist_l(n, 2), nlist_l(n, 4), nlist_r(n, 1), fx)
  fv(1:3) = fx * nlist_r(n, 2:4)
  forces(nlist_l(n, 1), :) = forces(nlist_l(n, 1), :) + fv(1:3)
  forces(nlist_l(n, 3), :) = forces(nlist_l(n, 3), :) - fv(1:3)
  
  
  ! A DENSITY AT B
  f = 1
  DO WHILE(fgroups_dens(nlist_l(n, 2), f) .NE. 0)
    fn = fgroups_dens(nlist_l(n, 2), f)
    CALL pot_search_b(2, nlist_l(n, 2), fn, nlist_r(n, 1), fx_d)  
    density(nlist_l(n, 3), fn) = density(nlist_l(n, 3), fn) + fx_d(1) 
    density_grad_ab(n, fn) = density_grad_ab(n, fn) + fx_d(2) 
    f = f + 1
  END DO 
  
  ! B DENSITY AT A
  f = 1
  DO WHILE(fgroups_dens(nlist_l(n, 4), f) .NE. 0)
    fn = fgroups_dens(nlist_l(n, 4), f)
    CALL pot_search_b(2, nlist_l(n, 4), fn, nlist_r(n, 1), fx_d)
    density(nlist_l(n, 1), fn) = density(nlist_l(n, 1), fn) + fx_d(1)
    density_grad_ba(n, fn) = density_grad_ba(n, fn) + fx_d(2) 
    f = f + 1
  END DO
END DO

DO n = 1, nl_count
  f = 1
  DO WHILE(fgroups_embe(nlist_l(n, 4), f) .NE. 0)
    fn = fgroups_embe(nlist_l(n, 2), f)
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_c(3, nlist_l(n, 2), fn, density(nlist_l(n, 1), fn), fx)
    epA = fx
  
    ! @Fi(p)/@p    Gradient of embedding energy for atom A
    CALL pot_search_c(3, nlist_l(n, 4), fn, density(nlist_l(n, 3), fn), fx)
    epB = fx
    
    ! @Pij(r)/@r   Gradient of density function at A due to atom B
    dpAB = density_grad_ab(n,fn)
      
    ! @Pji(r)/@r   Gradient of density function at B due to atom A
    dpBA = density_grad_ba(n,fn)

  
    fv(:) = (epA * dpBA + epB * dpAB) * nlist_r(n, 2:4)
    
    forces(nlist_l(n, 1), :) = forces(nlist_l(n, 1), :) + fv(1:3)
    forces(nlist_l(n, 3), :) = forces(nlist_l(n, 3), :) - fv(1:3)

    f = f + 1    
  END DO
END DO

DEALLOCATE(density)
DEALLOCATE(density_grad_ab)
DEALLOCATE(density_grad_ba)

END SUBROUTINE force_calc
