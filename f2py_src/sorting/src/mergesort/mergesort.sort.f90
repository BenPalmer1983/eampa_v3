SUBROUTINE ms_sort()
!#############################################################
! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
!#############################################################
IF(choice .EQ. 1)THEN
  CALL sort_1d_int()
ELSE IF(choice .EQ. 2) THEN
  CALL sort_2d_int()
ELSE IF(choice .EQ. 3) THEN
  CALL sort_1d_dp()
ELSE IF(choice .EQ. 4) THEN
  CALL sort_2d_dp()
END IF
END SUBROUTINE ms_sort
