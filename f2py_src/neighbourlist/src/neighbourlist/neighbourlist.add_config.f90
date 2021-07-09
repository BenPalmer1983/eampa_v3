SUBROUTINE add_config(rcut_in, a0_in, uv_in, copies_in, labels_in, coords_in)
!###########################################################
REAL(kind=DoubleReal) :: rcut_in
REAL(kind=DoubleReal) :: a0_in
REAL(kind=DoubleReal) :: uv_in(1:3,1:3)
INTEGER(kind=StandardInteger) :: copies_in(1:3)
INTEGER(kind=StandardInteger) :: labels_in(:)
REAL(kind=DoubleReal) :: coords_in(:,:)
!###########################################################
INTEGER(kind=StandardInteger) :: n, m, nx, ny, nz, k
REAL(kind=DoubleReal) :: copies(1:3,1:3)
!###########################################################

cc = cc + 1

copies = 0.0D0
copies(1,1) = copies_in(1)
copies(2,2) = copies_in(2)
copies(3,3) = copies_in(3)

rcut(cc) = rcut_in
a0(cc) = a0_in
uv(cc, 1:3, 1:3) = matmul(copies, uv_in)

IF(cc .EQ. 1)THEN
  n = 1
ELSE
  n = coords_key(cc - 1, 2) + 1
END IF

k = 0
coords_key(cc, 1) = n
DO nx = 1, copies_in(1)
  DO ny = 1, copies_in(2)
    DO nz = 1, copies_in(3)
      DO m = 1, SIZE(coords_in, 1)
        k = k + 1
        atom_id(n) = k
        atom_type(n) = labels_in(m)
        coords_crystal(n, 1) = (1.0D0 * (nx - 1) + coords_in(m, 1)) / (1.0D0 * copies_in(1))
	      coords_crystal(n, 2) = (1.0D0 * (ny - 1) + coords_in(m, 2)) / (1.0D0 * copies_in(2))
		    coords_crystal(n, 3) = (1.0D0 * (nz - 1) + coords_in(m, 3)) / (1.0D0 * copies_in(3))
        n = n + 1
      END DO
    END DO
  END DO
END DO
coords_key(cc, 2) = n - 1
coords_key(cc, 3) = n - coords_key(cc, 1) 


!###########################################################
END SUBROUTINE add_config