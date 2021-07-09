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

MODULE neighbourlist

USE kinds

IMPLICIT NONE

INCLUDE "neighbourlist.globals.f90"
INCLUDE "neighbourlist.interfaces.f90"

CONTAINS

INCLUDE "neighbourlist.init.f90"
INCLUDE "neighbourlist.add_config.f90"
INCLUDE "neighbourlist.build.f90"
INCLUDE "neighbourlist.make_real.f90"
INCLUDE "neighbourlist.make_halo.f90"
INCLUDE "neighbourlist.make_nl.f90"





END MODULE neighbourlist
