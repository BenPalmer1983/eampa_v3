!#############################################################
!#  HOW TO USE:
!#
!#
!# interp.interp4(x, x_points, y_points, y)               # Estimates f(x) at x    4 data points
!# interp.interp4dydx(x, x_points, y_points, y)           # Estimates f'(x) at x   4 data points
!# interp.interp4dydx(x, x_points, y_points, y)           # Estimates f^n(x) at x  4 data points
!#
!# interp.interpn(x, x_points, y_points, y)               # Estimates f(x) at x    n data points
!# interp.interp4dydx(x, x_points, y_points, y)           # Estimates f'(x) at x   n data points
!# interp.interp4dydx(x, x_points, y_points, y)           # Estimates f^p(x) at x  n data points
!#
!# interp.interpolate(xi, x, y, n_interp_in, yi)          # Searches large dataset, takes smaller subset and interpolates
!# interp.trap(xi, x, y, yi)                              # Very simple, searches dataset, line between 2 cloest points 
!#                                                        # used e.g. cross section data, not very well behaved data
!#
!# interp.fill(x, y, arr_l, arr_w, n_interp_in, out_arr)  # Takes data x, y and changes size in new array
!#                                                        # columns of new array x, f(x), f'(x), f''(x)...and so on
!#
!# fill_zoor(x, y, arr_l, arr_w,                          # Same as fill, but sets min/max for x
!#           x_min_range, x_max_range,                    # fills with zero out side of this range
!#           n_interp_in, out_arr)
!#                                                        
!#
!#############################################################

MODULE interp

USE kinds

IMPLICIT NONE

INCLUDE "interp.globals.f90"
INCLUDE "interp.interfaces.f90"

CONTAINS

INCLUDE "interp.interp4.f90"
INCLUDE "interp.interp4dydx.f90"
INCLUDE "interp.interp4dydxn.f90"

INCLUDE "interp.interpn.f90"
INCLUDE "interp.interpndydx.f90"
INCLUDE "interp.interpndydxn.f90"


INCLUDE "interp.interpolate.f90"
INCLUDE "interp.trap.f90"

INCLUDE "interp.fill.f90"
INCLUDE "interp.fill_zoor.f90"

END MODULE interp
