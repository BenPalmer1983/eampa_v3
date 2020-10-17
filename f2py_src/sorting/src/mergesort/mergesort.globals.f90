!#############################################################
!
!    GLOBAL VARIABLES
!
!#############################################################


REAL(kind=DoubleReal) ::                             start_time
REAL(kind=DoubleReal) ::                             end_time
REAL(kind=DoubleReal) ::                             time

INTEGER(KIND=StandardInteger) ::                     choice = 1      ! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
INTEGER(KIND=StandardInteger), ALLOCATABLE ::        key_table(:)
INTEGER(KIND=StandardInteger), ALLOCATABLE ::        arr_int_1d(:)
INTEGER(KIND=StandardInteger), ALLOCATABLE ::        arr_int_2d(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::                arr_dp_1d(:)
REAL(kind=DoubleReal), ALLOCATABLE ::                arr_dp_2d(:,:)

LOGICAL ::                                           sort_ascending = .TRUE.
LOGICAL ::                                           make_key_table = .TRUE.

INTEGER(KIND=StandardInteger) ::                     col = 0






