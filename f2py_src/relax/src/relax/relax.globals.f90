!###################################################################################
!#
!#                      ENERGIES GLOBALS
!#
!###################################################################################

!# MAX 1000 configs


!# CONFIGS + NL
INTEGER(kind=StandardInteger) ::                        atom_count = 0
REAL(kind=DoubleReal) ::                                rcut = 0.0D0
REAL(kind=DoubleReal) ::                                alat = 0.0D0
REAL(kind=DoubleReal) ::                                uv(1:3,1:3) = 0.0D0
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           labels(:)     
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords(:, :)  
REAL(kind=DoubleReal) ::                                gthreshold = 1.1D0 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostn(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostlabels(:)   
REAL(kind=DoubleReal), ALLOCATABLE ::                   ghostcoords(:,:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   ghostshift(:,:) 
LOGICAL, ALLOCATABLE ::                                 ghosthalo(:)  
INTEGER(kind=StandardInteger) ::                        ghost_count = 0
INTEGER(kind=StandardInteger) ::                        nl_count = 0
REAL(kind=DoubleReal), ALLOCATABLE ::                   nlist_r(:,:)  
INTEGER(kind=StandardInteger), ALLOCATABLE ::           nlist_l(:,:)  
LOGICAL, ALLOCATABLE ::                                 nlisthalo(:)    
REAL(kind=DoubleReal) ::                                to_cartesian(1:3,1:3) = 0.0D0  
REAL(kind=DoubleReal) ::                                from_cartesian(1:3,1:3) = 0.0D0  



!# POTENTIALS
INTEGER(kind=StandardInteger), PARAMETER ::             max_zbl = 100
INTEGER(kind=StandardInteger) ::                        pc = 0
INTEGER(kind=StandardInteger) ::                        mapkeys = 0
INTEGER(kind=StandardInteger) ::                        map_pn = 0
INTEGER(kind=StandardInteger) ::                        r_size = 1001
INTEGER(kind=StandardInteger) ::                        r_width = 4
INTEGER(kind=StandardInteger) ::                        pkey_temp(1:1000, 1:6)
INTEGER(kind=StandardInteger) ::                        pkey(1:1000, 1:6)
REAL(kind=DoubleReal) ::                                pot_rcut(1:1000) 
REAL(kind=DoubleReal) ::                                pot_temp(1:50000, 1:4) 
REAL(kind=DoubleReal) ::                                pot(1:50000, 1:4) 
REAL(kind=DoubleReal) ::                                pot_min_max(1:1000, 1:2) = -1.0D0
INTEGER(kind=StandardInteger) ::                        potmap(1:3, 1:50, 1:50) = 0         ! TYPE,  A,  B/Group
INTEGER(kind=StandardInteger) ::                        fgroups_dens(1:50, 1:51) = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
INTEGER(kind=StandardInteger) ::                        fgroups_embe(1:50, 1:51) = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
INTEGER(kind=StandardInteger) ::                        group_max = 0
INTEGER(kind=StandardInteger) ::                        zbl_counter = 0
INTEGER(kind=StandardInteger) ::                        zbl_i(1:max_zbl, 1:8)           ! 1 id_1  2 id_2  3 spline_type
REAL(kind=DoubleReal) ::                                zbl_r(1:max_zbl, 1:8)           ! 1 z1  2 z2  3 ra  4 rb
LOGICAL ::                                              zbl_l(1:max_zbl) = .TRUE.       ! on/off


!# CALCULATED
REAL(kind=DoubleReal), ALLOCATABLE ::                   forces(:, :) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   forces_last(:, :) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   velocities(:, :)   
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords_last(:, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords_relaxed(:, :)



!# MD XYZ
INTEGER(kind=StandardInteger) ::                        md_steps = 0
REAL(kind=DoubleReal), ALLOCATABLE ::                   md_xyz(:,:,:)



!###################################################################################