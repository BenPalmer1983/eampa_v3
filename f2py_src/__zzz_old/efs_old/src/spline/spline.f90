MODULE spline

USE kinds
USE interp, ONLY: interpolate, interp4, interp4dydx, interp4dydxn, fill
USE sls, ONLY: sls_solve

IMPLICIT NONE

INCLUDE "spline.interfaces.f90"
INCLUDE "spline.globals.f90"

CONTAINS

!    SPLINES
!    1 - poly3        p(x) = a + bx + cx^2 + dx^3
!    2 - poly5        p(x) = a + bx + cx^2 + dx^3 + ex^4 + fx^5
!    3 - exp3         p(x) = exp(a + bx + cx^2 + dx^3)
!    4 - exp5         p(x) = exp(a + bx + cx^2 + dx^3 + ex^4 + fx^5)


INCLUDE "spline.spline_functions.f90"
INCLUDE "spline.spline_ab.f90"
INCLUDE "spline.spline_array.f90"


END MODULE spline
