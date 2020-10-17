INTEGER(KIND=StandardInteger), ALLOCATABLE ::        arr_int_1d(:)
INTEGER(KIND=StandardInteger), ALLOCATABLE ::        arr_int_2d(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::                arr_dp_1d(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                arr_dp_2d(:,:)
LOGICAL ::                                           sort_ascending = .TRUE.
INTEGER(KIND=StandardInteger) ::                     col = 0
INTEGER(KIND=StandardInteger) ::                     mb_breakpoint = 100

INTEGER(KIND=StandardInteger), ALLOCATABLE ::        keytable(:)