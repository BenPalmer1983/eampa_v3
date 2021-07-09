!###################################################################################
!#
!#                      EFS GLOBALS
!#
!###################################################################################

!# MAX 1000 configs
INTEGER(kind=StandardInteger), PARAMETER ::             max_configs = 1000
INTEGER(kind=StandardInteger), PARAMETER ::             max_zbl = 100
INTEGER(kind=StandardInteger), PARAMETER ::             max_potfunctions = 1000


INTEGER(kind=StandardInteger) ::                        key(1:max_configs, 1:32)   ! 1 2:          coords start, coords end,  
                                                                                   ! 3 4:          gcoords start, gcoords end, 
                                                                                   ! 5 6:          nl start, nl end
                                                                                   ! 11, 12, 13:   e on, f on, stress on (for rss calc)  
                                                                                   ! 20 atom_count
!# CONFIGS + NL
INTEGER(kind=StandardInteger) ::                        cc = 0
INTEGER(kind=StandardInteger) ::                        total_atoms = 0
INTEGER(kind=StandardInteger) ::                        l_nl_size = 0
REAL(kind=DoubleReal), ALLOCATABLE ::                   rcut(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   alat(:)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   volume(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   uv(:,:) 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           labels(:)     
REAL(kind=DoubleReal), ALLOCATABLE ::                   coords(:, :)  
REAL(kind=DoubleReal), ALLOCATABLE ::                   energies(:)   
REAL(kind=DoubleReal), ALLOCATABLE ::                   forces(:, :, :)    
REAL(kind=DoubleReal), ALLOCATABLE ::                   stresses(:, :, :) 
REAL(kind=DoubleReal) ::                                gthreshold = 1.1D0 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostids(:)   
INTEGER(kind=StandardInteger), ALLOCATABLE ::           ghostlabels(:)   
REAL(kind=DoubleReal), ALLOCATABLE ::                   ghostcoords(:,:)  
LOGICAL, ALLOCATABLE ::                                 ghosthalo(:)    
REAL(kind=DoubleReal), ALLOCATABLE ::                   nlist_r(:,:)  
INTEGER(kind=StandardInteger), ALLOCATABLE ::           nlist_l(:,:)  
LOGICAL, ALLOCATABLE ::                                 nlisthalo(:)      







!# CALCULATED
REAL(kind=DoubleReal) ::                                config_energy(1:max_configs, 1:6) = 0.0D0 ! PAIR, EMBED, TOTAL  then per atom 4-6
REAL(kind=DoubleReal), ALLOCATABLE ::                   config_forces(:, :, :) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   config_stresses(:, :, :) 
INTEGER(kind=StandardInteger), ALLOCATABLE ::           stresses_calculated(:) 
REAL(kind=DoubleReal) ::                                max_rho = 0.0D0




!# RSS                                                                                    !             1  2  3  4                 5  6  7  8
REAL(kind=DoubleReal) ::                                config_rss(1:max_configs, 1:10)   ! UNWEIGHTED: E, F, S, Total   WEIGHTED: E, F, S, Total   
REAL(kind=DoubleReal) ::                                rss_weights(1:4) = 1.0D0   


REAL(kind=DoubleReal) ::                                energy_rss = 0.0D0
REAL(kind=DoubleReal) ::                                force_rss = 0.0D0
REAL(kind=DoubleReal) ::                                stress_rss = 0.0D0
REAL(kind=DoubleReal) ::                                total_rss = 0.0D0
REAL(kind=DoubleReal) ::                                energy_rss_weighted = 0.0D0
REAL(kind=DoubleReal) ::                                force_rss_weighted = 0.0D0
REAL(kind=DoubleReal) ::                                stress_rss_weighted = 0.0D0
REAL(kind=DoubleReal) ::                                total_rss_weighted = 0.0D0


REAL(kind=DoubleReal), ALLOCATABLE ::                   residuals(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::                   residuals_long(:) 
INTEGER(kind=StandardInteger) ::                        residuals_size = 0
INTEGER(kind=StandardInteger) ::                        residuals_long_size = 0


REAL(kind=DoubleReal) ::                                max_density = -1.0D99   


!# TIMERS    
REAL(kind=DoubleReal) ::                                timer_temp = 0.0D0
REAL(kind=DoubleReal) ::                                nl_timer = 0.0D0 
REAL(kind=DoubleReal) ::                                nl_timer_sum = 0.0D0 
REAL(kind=DoubleReal) ::                                e_timer = 0.0D0 
REAL(kind=DoubleReal) ::                                e_timer_sum = 0.0D0 
REAL(kind=DoubleReal) ::                                ef_timer = 0.0D0 
REAL(kind=DoubleReal) ::                                ef_timer_sum = 0.0D0 
REAL(kind=DoubleReal) ::                                efs_timer = 0.0D0 
REAL(kind=DoubleReal) ::                                efs_timer_sum = 0.0D0 
REAL(kind=DoubleReal) ::                                rss_timer = 0.0D0 
REAL(kind=DoubleReal) ::                                rss_timer_sum = 0.0D0 





!###################################################################################