SUBROUTINE init()

! DEALLOCATE
IF(ALLOCATED(ids))THEN
  DEALLOCATE(ids)
END IF
IF(ALLOCATED(labels))THEN
  DEALLOCATE(labels)
END IF
IF(ALLOCATED(coords))THEN
  DEALLOCATE(coords)
END IF
IF(ALLOCATED(coords_last))THEN
  DEALLOCATE(coords_last)
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
IF(ALLOCATED(ghostshift))THEN
  DEALLOCATE(ghostshift)
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


IF(ALLOCATED(forces))THEN
  DEALLOCATE(forces)
END IF
IF(ALLOCATED(forces_last))THEN
  DEALLOCATE(forces_last)
END IF
IF(ALLOCATED(velocities))THEN
  DEALLOCATE(velocities)
END IF





! ALLOCATE  
ALLOCATE(ids(1:2000000))                   ! 16MB
ALLOCATE(labels(1:2000000))                ! 16MB
ALLOCATE(coords(1:2000000,1:6))            ! 32MB
ALLOCATE(coords_last(1:2000000,1:3))            ! 32MB

ALLOCATE(ghostn(1:2000000))              ! 16MB
ALLOCATE(ghostids(1:2000000))              ! 16MB
ALLOCATE(ghostlabels(1:2000000))           ! 16MB
ALLOCATE(ghostcoords(1:2000000, 1:3))      ! 48MB
ALLOCATE(ghostshift(1:2000000, 1:3))      ! 48MB
ALLOCATE(ghosthalo(1:2000000))             ! 16MB

ALLOCATE(nlist_l(1:2000000, 1:6))
ALLOCATE(nlist_r(1:2000000, 1:4))
ALLOCATE(nlisthalo(1:2000000))

ALLOCATE(forces(1:2000000, 1:3))
ALLOCATE(forces_last(1:2000000, 1:3))
ALLOCATE(velocities(1:2000000, 1:6))
ALLOCATE(coords_relaxed(1:2000000, 1:3))




rcut = 0.0D0
alat = 0.0D0
uv = 0.0D0
ids = 0
labels = 0
coords = 0.0D0
coords_last = 0.0D0
ghostids = 0
ghostlabels = 0
ghostcoords = 0.0D0
ghostshift = 0.0D0
ghosthalo = .TRUE.
nlist_l = 0.0D0
nlist_r = 0.0D0
nlisthalo = .TRUE.
forces = 0.0D0
forces_last = 0.0D0
velocities = 0.0D0


END SUBROUTINE init 