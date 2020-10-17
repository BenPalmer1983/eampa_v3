SUBROUTINE zero()
!###########################################################
! MAKE CARTESIAN 
INTEGER(KIND=StandardInteger) :: n, i
REAL(kind=DoubleReal) :: x, y, z
!###########################################################


x = MINVAL(coords(1:atom_count, 4))
y = MINVAL(coords(1:atom_count, 5))
z = MINVAL(coords(1:atom_count, 6))
!print *, x, y, z

DO n = 1, atom_count
  coords(n, 4) = coords(n, 4) - x
  coords(n, 5) = coords(n, 5) - y
  coords(n, 6) = coords(n, 6) - z
  coords(n, 1:3) = matmul(from_cartesian, coords(n, 4:6))
END DO

DO n = 1, atom_count
  DO i = 1,3
    CALL Modulus(coords(n, i), 1.0D0, coords(n, i))
  END DO
  coords(n, 4:6) = matmul(to_cartesian, coords(n, 1:3))
END DO
END SUBROUTINE zero