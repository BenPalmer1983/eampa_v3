!##########################################################
! Add config and increments config counter cc
! Makes ghost cell


SUBROUTINE add_config(rcut_in, alat_in, uv_in, labels_in, coords_in)
!###########################################################
REAL(kind=DoubleReal) :: rcut_in
REAL(kind=DoubleReal) :: alat_in
REAL(kind=DoubleReal) :: uv_in(1:3,1:3)
INTEGER(kind=StandardInteger) :: labels_in(:)
REAL(kind=DoubleReal) :: coords_in(:,:)

!###########################################################
INTEGER(kind=StandardInteger) :: ccz
REAL(kind=DoubleReal) :: to_cart(1:3, 1:3)

INTEGER(KIND=StandardInteger) :: cn, an, cx, cy, cz, n, m, i

REAL(kind=DoubleReal) :: shift(1:3)
REAL(kind=DoubleReal) :: coords_crystal(1:3)
REAL(kind=DoubleReal) :: coords_cartesian(1:3)
REAL(kind=DoubleReal) :: half(1:3)
REAL(kind=DoubleReal) :: c(1:3)
REAL(kind=DoubleReal) :: ct(1:3)
REAL(kind=DoubleReal) :: td(1:3)
REAL(kind=DoubleReal) :: rd(1:3) 
REAL(kind=DoubleReal) :: r(1:3) 
REAL(kind=DoubleReal) :: vol_cell(1:3,1:3)
!###########################################################

atom_count = SIZE(labels_in, 1)
rcut = rcut_in
alat = alat_in
uv(1:3,1:3) = uv_in(1:3,1:3) 
labels(1:atom_count) = labels_in(:)
coords(1:atom_count, 1:3) = coords_in(:, 1:3)


to_cartesian(1:3,1:3) = alat * uv(1:3 , 1:3)  
CALL inv(to_cartesian, from_cartesian)


! MAKE CARTESIAN 
m = 0
DO n = 1, atom_count
  m = m + 1
  ids(n) = m
  DO i = 1,3
    CALL Modulus(coords(n, i), 1.0D0, coords(n, i))
  END DO
  coords(n, 4:6) = matmul(to_cartesian, coords(n, 1:3))
END DO  


END SUBROUTINE add_config 

