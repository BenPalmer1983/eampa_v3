!###################################################################################
!#
!#                      EFS GLOBALS
!#
!###################################################################################

!# MAX 1000 configs
INTEGER(kind=StandardInteger), PARAMETER ::             max_configs = 1000
INTEGER(kind=StandardInteger), PARAMETER ::             max_zbl = 100


!# EDITABLE IN PYTHON
INTEGER(kind=StandardInteger) ::                        eos_size = 10
REAL(kind=DoubleReal) ::                                eos_strain = 0.05D0   ! Total strain
INTEGER(kind=StandardInteger) ::                        eos_points = 0        ! 2*eos_size + 1
REAL(kind=DoubleReal) ::                                ec_strain = 0.05D0    ! Total strain
INTEGER(kind=StandardInteger) ::                        ec_size = 10          ! 2*eos_size + 1
INTEGER(kind=StandardInteger) ::                        ec_points = 0         ! ec_size




INTEGER(kind=StandardInteger) ::                        cc = 0
INTEGER(kind=StandardInteger) ::                        key(1:max_configs, 1:32)


REAL(kind=DoubleReal), ALLOCATABLE ::                   rcut(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   alat(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   volume(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   uv(:,:) 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           labels(:)     
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords(:, :)  






REAL(kind=DoubleReal) ::                                gthreshold = 1.1D0 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostlabels(:)   
REAL(kind=DoubleReal), ALLOCATABLE ::                   ghostcoords(:,:)  
LOGICAL, ALLOCATABLE ::                                 ghosthalo(:)              



REAL(kind=DoubleReal), ALLOCATABLE ::                   nlist_r(:,:)  
INTEGER(kind=StandardInteger), ALLOCATABLE ::           nlist_l(:,:)  
LOGICAL, ALLOCATABLE ::                                 nlisthalo(:)      




!# POTENTIALS


INTEGER(kind=StandardInteger) ::                        pc = 0
INTEGER(kind=StandardInteger) ::                        mapkeys = 0
INTEGER(kind=StandardInteger) ::                        map_pn = 0
INTEGER(kind=StandardInteger) ::                        r_size = 1001
INTEGER(kind=StandardInteger) ::                        r_width = 4
INTEGER(kind=StandardInteger) ::                        pkey_temp(1:max_configs, 1:6)
INTEGER(kind=StandardInteger) ::                        pkey(1:max_configs, 1:6)

REAL(kind=DoubleReal) ::                                pot_rcut(1:max_configs) 
REAL(kind=DoubleReal) ::                                pot_temp(1:50000, 1:4) 
REAL(kind=DoubleReal) ::                                pot(1:50000, 1:4) 
REAL(kind=DoubleReal) ::                                pot_min_max(1:max_configs, 1:2) = -1.0D0

INTEGER(kind=StandardInteger) ::                        potmap(1:3, 1:50, 1:50) = 0         ! TYPE,  A,  B/Group
INTEGER(kind=StandardInteger) ::                        fgroups_dens(1:50, 1:51) = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
INTEGER(kind=StandardInteger) ::                        fgroups_embe(1:50, 1:51) = 0        ! LIST OF FGROUPS FOR EACH TYPE  COUNTER IN 51
INTEGER(kind=StandardInteger) ::                        group_max = 0

INTEGER(kind=StandardInteger) ::                        zbl_counter = 0
INTEGER(kind=StandardInteger) ::                        zbl_i(1:max_zbl, 1:8)           ! 1 id_1  2 id_2  3 spline_type
REAL(kind=DoubleReal) ::                                zbl_r(1:max_zbl, 1:8)           ! 1 z1  2 z2  3 ra  4 rb
LOGICAL ::                                              zbl_l(1:max_zbl) = .TRUE.       ! on/off

!# CALCULATED
REAL(kind=DoubleReal) ::                                config_energy(1:max_configs, 1:6) = 0.0D0 ! PAIR, EMBED, TOTAL 


!# BP DATA
INTEGER(kind=StandardInteger) ::                        bp_configs_count = 0
INTEGER(kind=StandardInteger) ::                        bp_keys_i(1:100, 1:2) = -1
REAL(kind=DoubleReal) ::                                bp_keys_r(1:100, 1:2) = 0.0D0
INTEGER(kind=StandardInteger) ::                        bp_configs(1:100, 1:100) = 0

INTEGER(kind=StandardInteger) ::                        bp_config_data(1:100, 1:20) = 0     ! 1 eos_start, 2 eos_end, 3 ec1_start, 4 ec1_end...

REAL(kind=DoubleReal), ALLOCATABLE ::                   bp_results(:, :)    ! 1 v0 2 e0 3 b0 4 bop


REAL(kind=DoubleReal) ::                                bp_data_r(1:1000, 1:50) = 0.0D0     ! Strain etc  (doubles)
INTEGER(kind=StandardInteger) ::                        bp_data_i(1:1000, 1:50) = 0         ! Strain etc  (integers)


! INPUT 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           known_atoms_per_crystal(:)
INTEGER(kind=StandardInteger), ALLOCATABLE ::           known_expansion(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_alat(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_e0(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_b0(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_ec(:, :, :)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_g(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_e(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_v(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   known_amu_per_crystal(:)
LOGICAL, ALLOCATABLE ::                                 known_set(:,:)     ! 1 alat, 2 e0, 3 b0, 4 ec, 5 g, 6 e, 7 v


! CONFIG KEYS
INTEGER(kind=StandardInteger), ALLOCATABLE ::           bp_cn(:, :, :)

! FITTING ETC
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_strains(:, :, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_volumes(:, :, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_energies(:, :, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_energies_fit(:, :, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_fitting(:, :, :)
INTEGER(kind=StandardInteger), ALLOCATABLE ::           calc_sizes(:, :)

! RESULTS
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_alat(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_v0(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_e0(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_b0(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_rho(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_ec(:, :, :)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_sc(:, :, :) 



REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_b0_gpa(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_ec_gpa(:, :, :)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_sc_gpa(:, :, :) 

! BULK MODULUS
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_b0_r(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_b0_v(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_b0_avg(:)  

! SHEAR
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_g_r(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_g_v(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_g_avg(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_g_vec(:,:)  

! YOUNG
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_e(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_e_vec(:,:)  


! POISSON RATIO
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_v(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_v_vec(:, :, :)  


! Velocities
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_vl(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_vt(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_vm(:) 

! Temperatures
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_melting(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   calc_debye(:) 



REAL(kind=DoubleReal) ::                                rss_weights(1:10) = 1.0D0  

REAL(kind=DoubleReal), ALLOCATABLE ::                   rss(:, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   rss_total(:)
REAL(kind=DoubleReal) ::                                rss_total_rss = 0.0D0

REAL(kind=DoubleReal), ALLOCATABLE ::                   rss_w(:, :)
REAL(kind=DoubleReal), ALLOCATABLE ::                   rss_total_w(:)
REAL(kind=DoubleReal) ::                                rss_total_rss_w = 0.0D0



!###################################################################################