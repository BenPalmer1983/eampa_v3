SUBROUTINE run_md(steps, dt)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: steps
REAL(kind=DoubleReal), INTENT(IN) :: dt
!###########################################################
INTEGER(KIND=StandardInteger) :: n, step
REAL(kind=DoubleReal) :: limit = 1.0E-4
REAL(kind=DoubleReal) :: v(1:3)
!###########################################################

IF(ALLOCATED(md_xyz))THEN
  DEALLOCATE(md_xyz)
END IF
ALLOCATE(md_xyz(1:steps+1, 1:atom_count, 1:3))

velocities(1:atom_count, :) = 0.0D0

 
md_xyz(1, 1:atom_count, 1:3) = coords(1:atom_count, 4:6)
md_steps = steps

CALL make_ghost() 
CALL make_nl()             
CALL force_calc()

DO step =1, steps
  ! Calc V
  velocities(1:atom_count, 4:6) = velocities(1:atom_count, 1:3) + 0.5D0 * dt * forces(1:atom_count, 1:3)
  
  ! Calc x new
  coords(1:atom_count, 4:6) = coords(1:atom_count, 4:6) + dt * velocities(1:atom_count, 4:6)
  md_xyz(step + 1, 1:atom_count, 1:3) = coords(1:atom_count, 4:6)
  
  CALL refresh_coords()
  CALL ghost_update() 
  CALL nl_update()             
  CALL force_calc()
  
  velocities(1:atom_count, 1:3) = velocities(1:atom_count, 4:6) + 0.5 * dt * forces(1:atom_count, 1:3)
  !print *, coords(1, 4:6)

END DO




!###########################################################
END SUBROUTINE run_md















SUBROUTINE run_md_old(steps, dt)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: steps
REAL(kind=DoubleReal), INTENT(IN) :: dt
!###########################################################
INTEGER(KIND=StandardInteger) :: n, step
REAL(kind=DoubleReal) :: limit = 1.0E-4
REAL(kind=DoubleReal) :: v(1:3)
!###########################################################

IF(ALLOCATED(md_xyz))THEN
  DEALLOCATE(md_xyz)
END IF
ALLOCATE(md_xyz(1:steps+1, 1:atom_count, 1:3))

velocities(1:atom_count, :) = 0.0D0

CALL make_ghost() 
CALL make_nl()             
CALL force_calc()

DO step =1, steps
  ! Calc V
  velocities(1:atom_count, 4:6) = velocities(1:atom_count, 1:3) + 0.5D0 * dt * forces(1:atom_count, 1:3)
  
  ! Calc x new
  coords(1:atom_count, 4:6) = coords(1:atom_count, 4:6) + dt * velocities(1:atom_count, 4:6)
  md_xyz(step + 1, 1:atom_count, 1:3) = coords(1:atom_count, 4:6)
  
  CALL refresh_coords()
  CALL make_ghost() 
  CALL make_nl()             
  CALL force_calc()
  
  velocities(1:atom_count, 1:3) = velocities(1:atom_count, 4:6) + 0.5 * dt * forces(1:atom_count, 1:3)
  print *, coords(1, 4:6)

END DO




!###########################################################
END SUBROUTINE run_md_old