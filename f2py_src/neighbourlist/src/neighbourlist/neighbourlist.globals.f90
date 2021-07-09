

INTEGER(kind=StandardInteger), PARAMETER ::             max_configs = 1000
INTEGER(kind=StandardInteger), PARAMETER ::             max_coords =  1000000
INTEGER(kind=StandardInteger), PARAMETER ::             max_coords_halo =  4000000
INTEGER(kind=StandardInteger), PARAMETER ::             max_nl =  10000000



INTEGER(kind=StandardInteger) ::                        cc = 0


REAL(kind=DoubleReal) ::                                rcut(1:max_configs)
REAL(kind=DoubleReal) ::                                a0(1:max_configs)
REAL(kind=DoubleReal) ::                                uv(1:max_configs,1:3,1:3)

! KEYS
INTEGER(kind=StandardInteger) ::                        coords_key(1:max_configs, 1:3)
INTEGER(kind=StandardInteger) ::                        halo_key(1:max_configs, 1:3)
INTEGER(kind=StandardInteger) ::                        nl_key(1:max_configs, 1:3)

! COORDS
INTEGER(kind=StandardInteger), ALLOCATABLE ::           atom_id(:) 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           atom_type(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords_crystal(:, :) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords_real(:, :) 

! HALO
INTEGER(kind=StandardInteger), ALLOCATABLE ::           atom_halo_id(:) 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           atom_halo_type(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords_halo(:, :)  
LOGICAL, ALLOCATABLE ::                                 in_halo(:) 

! NEIGHBOUR LIST
INTEGER(kind=StandardInteger), ALLOCATABLE ::           nl_atoms(:,:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   nl_r(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   nl_rvec(:, :)  
LOGICAL, ALLOCATABLE ::                                 nl_in_halo(:) 

