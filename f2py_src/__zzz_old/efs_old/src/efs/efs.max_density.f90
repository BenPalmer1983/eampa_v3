SUBROUTINE max_density_calc()
!###########################################################
INTEGER(kind=StandardInteger) :: pn
REAL(kind=DoubleReal) :: fx
REAL(kind=DoubleReal) :: temp_rho
INTEGER(kind=StandardInteger) :: i, j, k, n, m, c

REAL(kind=DoubleReal) :: r(1:100)
INTEGER(kind=StandardInteger) :: rn(1:100)
!###########################################################


!# FCC
r(1) = 7.48332D0
r(2) = 6.32456D0
r(3) = 6.92821D0
r(4) = 4.89898D0
r(5) = 5.65686D0
r(6) = 2.82843D0
r(7) = 4.0D0
rn(1) = 48
rn(2) = 24
rn(3) = 8
rn(4) = 24
rn(5) = 12
rn(6) = 12
rn(7) = 6

temp_rho = 0.0D0
DO m=1, 7
  CALL pot_search_a(2, 1, 1, r(m), fx)
  temp_rho = temp_rho + (1.0D0 * rn(m)) * fx
END DO
print *, temp_rho



END SUBROUTINE max_density_calc
