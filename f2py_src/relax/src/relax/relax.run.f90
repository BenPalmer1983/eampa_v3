SUBROUTINE run(steps, dt, damp)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: steps
REAL(kind=DoubleReal), INTENT(IN) :: dt
REAL(kind=DoubleReal), INTENT(IN) :: damp
!###########################################################
INTEGER(KIND=StandardInteger) :: n, step, i, dn
REAL(kind=DoubleReal) :: limit = 1.0E-4
REAL(kind=DoubleReal) :: v(1:3)
REAL(kind=DoubleReal) :: alpha = 0.0D0
REAL(kind=DoubleReal) :: frss = 0.0D0
REAL(kind=DoubleReal) :: frss_last = 0.0D0
!###########################################################


velocities(1:atom_count, :) = 0.0D0

CALL make_ghost() 
CALL make_nl()             
CALL force_calc()
frss = sum(abs(forces(1:atom_count, 1:3)))
print *, frss

alpha = 0.1D0
DO step =1, 1   
  coords_last(1:atom_count, 1:3) = coords(1:atom_count, 4:6)
  forces_last(1:atom_count, 1:3) = forces(1:atom_count, 1:3) 
  
  DO dn=1,5
    coords(1:atom_count, 4:6) = alpha * forces_last(1:atom_count, 1:3) +  coords_last(1:atom_count, 1:3)
      
    CALL refresh_coords()
    CALL ghost_update() 
    CALL nl_update()  
    CALL force_calc()
    frss = sum(abs(forces(1:atom_count, 1:3)))
    alpha = 0.5D0 * alpha
    IF(dn .GT. 1 .AND. frss .GT. frss_last)THEN
      print *, frss_last, frss, "BREAK"
    END IF
    frss_last = frss
  END DO
END DO


CALL refresh_coords()
CALL ghost_update() 
CALL nl_update()  
CALL force_calc()
frss = sum(abs(forces(1:atom_count, 1:3)))
print *, frss


CALL zero()

!###########################################################
END SUBROUTINE run


SUBROUTINE run_damped(steps, dt, damp)
!###########################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: steps
REAL(kind=DoubleReal), INTENT(IN) :: dt
REAL(kind=DoubleReal), INTENT(IN) :: damp
!###########################################################
INTEGER(KIND=StandardInteger) :: n, step, i, dn
REAL(kind=DoubleReal) :: limit = 1.0E-4
REAL(kind=DoubleReal) :: v(1:3)
REAL(kind=DoubleReal) :: frss = 0.0D0
REAL(kind=DoubleReal) :: frss_last = 0.0D0
!###########################################################


velocities(1:atom_count, :) = 0.0D0

CALL make_ghost() 
CALL make_nl()             
CALL force_calc()
frss = sum(abs(forces(1:atom_count, 1:3)))

DO step =1, steps 
  ! Store last force
  forces_last(1:atom_count, 1:3) = forces(1:atom_count, 1:3) 
  frss_last = frss
  
  ! Calc V
  velocities(1:atom_count, 4:6) = velocities(1:atom_count, 1:3) + 0.5D0 * dt * forces(1:atom_count, 1:3)
  
  ! Calc x new
  coords_last(1:atom_count, 1:3) = coords(1:atom_count, 4:6)
  coords(1:atom_count, 4:6) = coords(1:atom_count, 4:6) + dt * velocities(1:atom_count, 4:6)
  
  CALL refresh_coords()
  CALL ghost_update() 
  CALL nl_update()  
  CALL force_calc()
  frss = sum(abs(forces(1:atom_count, 1:3)))
  
  velocities(1:atom_count, 1:3) = velocities(1:atom_count, 4:6) + 0.5 * dt * forces(1:atom_count, 1:3)

  DO n=1,atom_count
    DO i=1,3
      IF(forces_last(n, i) * forces(n, i) .LT. 0.0D0) THEN
        velocities(n, i) = damp * velocities(n, i)
      END IF
    END DO
  END DO
END DO

CALL refresh_coords()
CALL ghost_update() 
CALL nl_update()  
CALL force_calc()
frss = sum(abs(forces(1:atom_count, 1:3)))

CALL zero()

!###########################################################
END SUBROUTINE run_damped