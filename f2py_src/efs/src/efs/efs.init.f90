SUBROUTINE init()

!# RESET KEY + CONFIG COUNTER
key = 0
cc = 0


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
IF(ALLOCATED(energies))THEN
  DEALLOCATE(energies)
END IF
IF(ALLOCATED(forces))THEN
  DEALLOCATE(forces)
END IF
IF(ALLOCATED(stresses))THEN
  DEALLOCATE(stresses)
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

IF(ALLOCATED(config_forces))THEN
  DEALLOCATE(config_forces)
END IF
IF(ALLOCATED(config_stresses))THEN
  DEALLOCATE(config_stresses)
END IF



! ALLOCATE  
ALLOCATE(rcut(1:1000))                     ! 8KB
ALLOCATE(alat(1:1000))                     ! 8KB
ALLOCATE(volume(1:1000))                   ! 8KB
ALLOCATE(uv(1:3000, 1:3))                  ! 74KB
ALLOCATE(ids(1:2000000))                   ! 16MB
ALLOCATE(labels(1:2000000))                ! 16MB
ALLOCATE(coords(1:2000000,1:6))            ! 32MB
ALLOCATE(energies(1:2000000))             
ALLOCATE(forces(1:2000000,1:3))            ! 16MB
ALLOCATE(stresses(1:3000,1:3))             ! 74KB

ALLOCATE(ghostids(1:2000000))              ! 16MB
ALLOCATE(ghostlabels(1:2000000))           ! 16MB
ALLOCATE(ghostcoords(1:2000000, 1:3))      ! 48MB
ALLOCATE(ghosthalo(1:2000000))             ! 16MB

ALLOCATE(nlist_l(1:2000000, 1:4))
ALLOCATE(nlist_r(1:2000000, 1:4))
ALLOCATE(nlisthalo(1:2000000))

ALLOCATE(config_forces(1:2000000,1:3))     ! 16MB
ALLOCATE(config_stresses(1:3000,1:3))             ! 74KB

rcut = 0.0D0
alat = 0.0D0
volume = 0.0D0
uv = 0.0D0
ids = 0
labels = 0
coords = 0.0D0
energies = 0.0D0
forces = 0.0D0
stresses = 0.0D0
ghostids = 0
ghostlabels = 0
ghostcoords = 0.0D0
ghosthalo = .TRUE.
nlist_l = 0.0D0
nlist_r = 0.0D0
nlisthalo = .TRUE.
config_forces = 0.0D0
config_stresses = 0.0D0





END SUBROUTINE init 