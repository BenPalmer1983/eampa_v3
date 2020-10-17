PROGRAM chi2prog
! University of Birmingham
! Ben Palmer
! p(chisq) function from https://www.mathsisfun.com/data/images/chi-square.js

USE kinds
USE chi2, ONLY: run

IMPLICIT NONE

CALL main()

CONTAINS

SUBROUTINE main()
!#################################
CHARACTER(LEN=128) :: filename
CALL get_command_argument(1, filename)
CALL run(TRIM(ADJUSTL(filename)))
END SUBROUTINE main


END PROGRAM chi2prog
