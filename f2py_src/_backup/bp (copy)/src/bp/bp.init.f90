SUBROUTINE init(v_config_size_in, v_ghost_size_in, v_nl_size_in)
!###########################################################
INTEGER(KIND=StandardInteger), OPTIONAL :: v_config_size_in
INTEGER(KIND=StandardInteger), OPTIONAL :: v_ghost_size_in
INTEGER(KIND=StandardInteger), OPTIONAL :: v_nl_size_in
!###########################################################

IF(PRESENT(v_config_size_in))THEN
  v_config_size = v_config_size_in
END IF
IF(PRESENT(v_ghost_size_in))THEN
  v_ghost_size = v_ghost_size_in
END IF
IF(PRESENT(v_nl_size_in))THEN
  v_nl_size = v_nl_size_in
END IF



!# RESET KEY + CONFIG COUNTER
key = 0
cc = 0
cc_log = 0






! DEALLOCATE
IF(ALLOCATED(rcut))THEN
  DEALLOCATE(rcut)
END IF
IF(ALLOCATED(alat))THEN
  DEALLOCATE(alat)
END IF
IF(ALLOCATED(volume))THEN
  DEALLOCATE(volume)
END IF
IF(ALLOCATED(uv))THEN
  DEALLOCATE(uv)
END IF
IF(ALLOCATED(labels))THEN
  DEALLOCATE(labels)
END IF
IF(ALLOCATED(coords))THEN
  DEALLOCATE(coords)
END IF


IF(ALLOCATED(ghostids))THEN
  DEALLOCATE(ghostids)
END IF
IF(ALLOCATED(ghostlabels))THEN
  DEALLOCATE(ghostlabels)
END IF
IF(ALLOCATED(ghostcoords))THEN
  DEALLOCATE(ghostcoords)
END IF
IF(ALLOCATED(ghosthalo))THEN
  DEALLOCATE(ghosthalo)
END IF

IF(ALLOCATED(nlist_l))THEN
  DEALLOCATE(nlist_l)
END IF
IF(ALLOCATED(nlist_r))THEN
  DEALLOCATE(nlist_r)
END IF
IF(ALLOCATED(nlisthalo))THEN
  DEALLOCATE(nlisthalo)
END IF


IF(ALLOCATED(known_atoms_per_crystal))THEN
  DEALLOCATE(known_atoms_per_crystal)
END IF
IF(ALLOCATED(known_expansion))THEN
  DEALLOCATE(known_expansion)
END IF
IF(ALLOCATED(known_amu_per_crystal))THEN
  DEALLOCATE(known_amu_per_crystal)
END IF
IF(ALLOCATED(known_alat))THEN
  DEALLOCATE(known_alat)
END IF
IF(ALLOCATED(known_b0))THEN
  DEALLOCATE(known_b0)
END IF
IF(ALLOCATED(known_e0))THEN
  DEALLOCATE(known_e0)
END IF
IF(ALLOCATED(known_ec))THEN
  DEALLOCATE(known_ec)
END IF
IF(ALLOCATED(known_g))THEN
  DEALLOCATE(known_g)
END IF
IF(ALLOCATED(known_e))THEN
  DEALLOCATE(known_e)
END IF
IF(ALLOCATED(known_v))THEN
  DEALLOCATE(known_v)
END IF


IF(ALLOCATED(known_set))THEN
  DEALLOCATE(known_set)
END IF


IF(ALLOCATED(bp_cn))THEN
  DEALLOCATE(bp_cn)
END IF
IF(ALLOCATED(bp_results))THEN
  DEALLOCATE(bp_results)
END IF


IF(ALLOCATED(calc_fitting))THEN
  DEALLOCATE(calc_fitting)
END IF
IF(ALLOCATED(calc_strains))THEN
  DEALLOCATE(calc_strains)
END IF
IF(ALLOCATED(calc_energies))THEN
  DEALLOCATE(calc_energies)
END IF
IF(ALLOCATED(calc_energies_fit))THEN
  DEALLOCATE(calc_energies_fit)
END IF
IF(ALLOCATED(calc_volumes))THEN
  DEALLOCATE(calc_volumes)
END IF
IF(ALLOCATED(calc_sizes))THEN
  DEALLOCATE(calc_sizes)
END IF

IF(ALLOCATED(calc_alat))THEN
  DEALLOCATE(calc_alat)
END IF
IF(ALLOCATED(calc_v0))THEN
  DEALLOCATE(calc_v0)
END IF
IF(ALLOCATED(calc_e0))THEN
  DEALLOCATE(calc_e0)
END IF
IF(ALLOCATED(calc_b0))THEN
  DEALLOCATE(calc_b0)
END IF
IF(ALLOCATED(calc_rho))THEN
  DEALLOCATE(calc_rho)
END IF
IF(ALLOCATED(calc_ec))THEN
  DEALLOCATE(calc_ec)
END IF
IF(ALLOCATED(calc_sc))THEN
  DEALLOCATE(calc_sc)
END IF
IF(ALLOCATED(calc_ec_gpa))THEN
  DEALLOCATE(calc_ec_gpa)
END IF
IF(ALLOCATED(calc_b0_gpa))THEN
  DEALLOCATE(calc_b0_gpa)
END IF
IF(ALLOCATED(calc_sc_gpa))THEN
  DEALLOCATE(calc_sc_gpa)
END IF


IF(ALLOCATED(calc_b0_r))THEN
  DEALLOCATE(calc_b0_r)
END IF
IF(ALLOCATED(calc_b0_v))THEN
  DEALLOCATE(calc_b0_v)
END IF
IF(ALLOCATED(calc_b0_v))THEN
  DEALLOCATE(calc_b0_avg)
END IF

IF(ALLOCATED(calc_g_r))THEN
  DEALLOCATE(calc_g_r)
END IF
IF(ALLOCATED(calc_g_v))THEN
  DEALLOCATE(calc_g_v)
END IF
IF(ALLOCATED(calc_g_avg))THEN
  DEALLOCATE(calc_g_avg)
END IF
IF(ALLOCATED(calc_g_vec))THEN
  DEALLOCATE(calc_g_vec)
END IF

IF(ALLOCATED(calc_e))THEN
  DEALLOCATE(calc_e)
END IF
IF(ALLOCATED(calc_e_vec))THEN
  DEALLOCATE(calc_e_vec)
END IF

IF(ALLOCATED(calc_v))THEN
  DEALLOCATE(calc_v)
END IF
IF(ALLOCATED(calc_v_vec))THEN
  DEALLOCATE(calc_v_vec)
END IF

IF(ALLOCATED(calc_vl))THEN
  DEALLOCATE(calc_vl)
END IF
IF(ALLOCATED(calc_vt))THEN
  DEALLOCATE(calc_vt)
END IF
IF(ALLOCATED(calc_vm))THEN
  DEALLOCATE(calc_vm)
END IF

IF(ALLOCATED(calc_melting))THEN
  DEALLOCATE(calc_melting)
END IF
IF(ALLOCATED(calc_debye))THEN
  DEALLOCATE(calc_debye)
END IF



IF(ALLOCATED(calc_b0_r_gpa))THEN
  DEALLOCATE(calc_b0_r_gpa)
END IF
IF(ALLOCATED(calc_b0_v_gpa))THEN
  DEALLOCATE(calc_b0_v_gpa)
END IF
IF(ALLOCATED(calc_g_r_gpa))THEN
  DEALLOCATE(calc_g_r_gpa)
END IF
IF(ALLOCATED(calc_g_v_gpa))THEN
  DEALLOCATE(calc_g_v_gpa)
END IF
IF(ALLOCATED(calc_g_avg_gpa))THEN
  DEALLOCATE(calc_g_avg_gpa)
END IF
IF(ALLOCATED(calc_g_avg_gpa))THEN
  DEALLOCATE(calc_g_avg_gpa)
END IF
IF(ALLOCATED(calc_g_vec_gpa))THEN
  DEALLOCATE(calc_g_vec_gpa)
END IF
IF(ALLOCATED(calc_e_gpa))THEN
  DEALLOCATE(calc_e_gpa)
END IF
IF(ALLOCATED(calc_e_vec_gpa))THEN
  DEALLOCATE(calc_e_vec_gpa)
END IF








IF(ALLOCATED(rss))THEN
  DEALLOCATE(rss)
END IF
IF(ALLOCATED(rss_total))THEN
  DEALLOCATE(rss_total)
END IF

IF(ALLOCATED(rss_w))THEN
  DEALLOCATE(rss_w)
END IF
IF(ALLOCATED(rss_total_w))THEN
  DEALLOCATE(rss_total_w)
END IF

IF(ALLOCATED(residuals))THEN
  DEALLOCATE(residuals)
END IF



! ALLOCATE  
ALLOCATE(rcut(1:1000))                     ! 8KB
ALLOCATE(alat(1:1000))                     ! 8KB
ALLOCATE(volume(1:1000))                     ! 8KB
ALLOCATE(uv(1:3000, 1:3))                  ! 74KB

! COORDS
ALLOCATE(ids(1:v_config_size))                   !
ALLOCATE(labels(1:v_config_size))                !
ALLOCATE(coords(1:v_config_size,1:6))            ! 

! GHOST COORDS
ALLOCATE(ghostids(1:v_ghost_size))              !
ALLOCATE(ghostlabels(1:v_ghost_size))           !
ALLOCATE(ghostcoords(1:v_ghost_size, 1:3))      !v_nl_size_in
ALLOCATE(ghosthalo(1:v_ghost_size))             !

ALLOCATE(nlist_l(1:v_nl_size, 1:4))
ALLOCATE(nlist_r(1:v_nl_size, 1:4))
ALLOCATE(nlisthalo(1:v_nl_size))


ALLOCATE(known_atoms_per_crystal(1:max_bps))
ALLOCATE(known_expansion(1:max_bps))
ALLOCATE(known_amu_per_crystal(1:max_bps))
ALLOCATE(known_alat(1:max_bps))
ALLOCATE(known_b0(1:max_bps))
ALLOCATE(known_e0(1:max_bps))
ALLOCATE(known_ec(1:max_bps, 1:6, 1:6))
ALLOCATE(known_g(1:max_bps))
ALLOCATE(known_e(1:max_bps))
ALLOCATE(known_v(1:max_bps))
ALLOCATE(known_set(1:max_bps, 1:20))

ALLOCATE(bp_cn(1:max_bps, 1:10, 1:2))
ALLOCATE(bp_results(1:max_bps, 1:100))



ALLOCATE(calc_fitting(1:max_bps, 1:9, 1:3))
ALLOCATE(calc_strains(1:max_bps, 1:10, 1:100))
ALLOCATE(calc_volumes(1:max_bps, 1:10, 1:100))
ALLOCATE(calc_energies(1:max_bps, 1:10, 1:100))
ALLOCATE(calc_energies_fit(1:max_bps, 1:10, 1:100))
ALLOCATE(calc_sizes(1:max_bps, 1:10))
ALLOCATE(calc_alat(1:max_bps))
ALLOCATE(calc_v0(1:max_bps))
ALLOCATE(calc_e0(1:max_bps))
ALLOCATE(calc_b0(1:max_bps))
ALLOCATE(calc_rho(1:max_bps))
ALLOCATE(calc_ec(1:max_bps, 1:6, 1:6))
ALLOCATE(calc_sc(1:max_bps, 1:6, 1:6))
ALLOCATE(calc_b0_gpa(1:max_bps))
ALLOCATE(calc_ec_gpa(1:max_bps, 1:6, 1:6))
ALLOCATE(calc_sc_gpa(1:max_bps, 1:6, 1:6))

ALLOCATE(calc_b0_r(1:max_bps))
ALLOCATE(calc_b0_v(1:max_bps))
ALLOCATE(calc_b0_avg(1:max_bps))

ALLOCATE(calc_g_r(1:max_bps))
ALLOCATE(calc_g_v(1:max_bps))
ALLOCATE(calc_g_avg(1:max_bps))
ALLOCATE(calc_g_vec(1:max_bps,1:3))

ALLOCATE(calc_e(1:max_bps))
ALLOCATE(calc_e_vec(1:max_bps,1:3))

ALLOCATE(calc_v(1:max_bps))
ALLOCATE(calc_v_vec(1:max_bps,1:3,1:3))

ALLOCATE(calc_vl(1:max_bps))
ALLOCATE(calc_vt(1:max_bps))
ALLOCATE(calc_vm(1:max_bps))

ALLOCATE(calc_melting(1:max_bps))
ALLOCATE(calc_debye(1:max_bps))


ALLOCATE(calc_b0_r_gpa(1:max_bps))
ALLOCATE(calc_b0_v_gpa(1:max_bps))
ALLOCATE(calc_b0_avg_gpa(1:max_bps))

ALLOCATE(calc_g_r_gpa(1:max_bps))
ALLOCATE(calc_g_v_gpa(1:max_bps))
ALLOCATE(calc_g_avg_gpa(1:max_bps))
ALLOCATE(calc_g_vec_gpa(1:max_bps,1:3))

ALLOCATE(calc_e_gpa(1:max_bps))
ALLOCATE(calc_e_vec_gpa(1:max_bps,1:3))



ALLOCATE(rss(1:max_bps, 1:20))
ALLOCATE(rss_total(1:max_bps))

ALLOCATE(rss_w(1:max_bps, 1:20))
ALLOCATE(rss_total_w(1:max_bps))



ALLOCATE(residuals(1:max_bps * 40))



rcut = 0.0D0
alat = 0.0D0
volume = 0.0D0
uv = 0.0D0
ids = 0
labels = 0
coords = 0.0D0
ghostids = 0
ghostlabels = 0
ghostcoords = 0.0D0
ghosthalo = .TRUE.
nlist_l = 0.0D0
nlist_r = 0.0D0
nlisthalo = .TRUE.
known_atoms_per_crystal = 0
known_expansion = 0
known_alat = 0.0D0
known_b0 = 0.0D0
known_e0 = 0.0D0
known_e0 = 0.0D0
known_g = 0.0D0
known_e = 0.0D0
known_v = 0.0D0
known_amu_per_crystal = 0.0D0
known_set = .FALSE.
bp_cn = 0
bp_results = 0.0D0
calc_fitting = 0.0D0
calc_strains = 0.0D0
calc_volumes = 0.0D0
calc_energies = 0.0D0
calc_energies_fit = 0.0D0
calc_sizes = 0
calc_alat = 0.0D0
calc_v0 = 0.0D0
calc_e0 = 0.0D0
calc_b0 = 0.0D0
calc_rho = 0.0D0
calc_ec = 0.0D0
calc_sc = 0.0D0
calc_b0_gpa = 0.0D0
calc_ec_gpa = 0.0D0
calc_sc_gpa = 0.0D0

calc_b0_r = 0.0D0
calc_b0_v = 0.0D0
calc_b0_avg = 0.0D0

calc_g_r = 0.0D0
calc_g_v = 0.0D0
calc_g_avg = 0.0D0
calc_g_vec = 0.0D0

calc_e = 0.0D0
calc_e_vec = 0.0D0

calc_v = 0.0D0
calc_v_vec = 0.0D0

calc_vl = 0.0D0
calc_vt = 0.0D0
calc_vm = 0.0D0

calc_melting = 0.0D0
calc_debye = 0.0D0

rss = 0.0D0
rss_total = 0.0D0

rss_w = 0.0D0
rss_total_w = 0.0D0


END SUBROUTINE init 