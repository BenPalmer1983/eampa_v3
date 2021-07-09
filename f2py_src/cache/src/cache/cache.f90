!#############################################################
!#  HOW TO USE:
!#
!#
!# 
!# 
!#
!#
!#
!#
!#############################################################

MODULE cache

USE kinds

IMPLICIT NONE

INCLUDE "cache.globals.f90"
INCLUDE "cache.interfaces.f90"

CONTAINS

INCLUDE "cache.init.f90"
INCLUDE "cache.clear.f90"
INCLUDE "cache.set.f90"
INCLUDE "cache.get.f90"




END MODULE cache
