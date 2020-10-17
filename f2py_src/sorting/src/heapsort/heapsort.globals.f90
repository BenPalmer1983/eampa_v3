!#############################################################
!
!    GLOBAL VARIABLES
!
!#############################################################

INTEGER(KIND=StandardInteger) ::               rows = 0
INTEGER(KIND=StandardInteger) ::               cols = 0

INTEGER(KIND=StandardInteger) ::               degrees_of_freedom = 0

REAL(kind=DoubleReal), ALLOCATABLE ::          observed(:,:) 
REAL(kind=DoubleReal), ALLOCATABLE ::          expected(:,:)

REAL(kind=DoubleReal), ALLOCATABLE ::          row_sums(:) 
REAL(kind=DoubleReal), ALLOCATABLE ::          col_sums(:) 
REAL(kind=DoubleReal) ::                       total

REAL(kind=DoubleReal), ALLOCATABLE ::          o_e(:)

REAL(kind=DoubleReal) ::                       chi_squared = 0.0D0
REAL(kind=DoubleReal) ::                       p_chi_squared = 0.0D0

REAL(kind=DoubleReal) ::                       log_sqrt_pi = 0.5723649429247000870717135D0
REAL(kind=DoubleReal) ::                       i_sqrt_pi = 0.5641895835477562869480795D0

REAL(kind=DoubleReal) ::                       big_x = 20.0D0