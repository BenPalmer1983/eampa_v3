SUBROUTINE ascending()
sort_ascending = .TRUE.
END SUBROUTINE ascending


SUBROUTINE descending()
sort_ascending = .TRUE.
END SUBROUTINE descending


SUBROUTINE set_1d_int(arr_in)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: arr_in(:)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
choice = 1      ! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
IF(ALLOCATED(arr_int_1d))THEN
  DEALLOCATE(arr_int_1d)
END IF
IF(ALLOCATED(key_table))THEN
  DEALLOCATE(key_table)
END IF
ALLOCATE(arr_int_1d(1: SIZE(arr_in,1)))
ALLOCATE(key_table(1: SIZE(arr_in,1)))
arr_int_1d(:) = arr_in(:)
DO n = 1, SIZE(arr_in,1)
  key_table(n) = n
END DO
END SUBROUTINE set_1d_int


SUBROUTINE set_2d_int(arr_in)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: arr_in(:,:)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
choice = 2      ! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
IF(ALLOCATED(arr_int_2d))THEN
  DEALLOCATE(arr_int_2d)
END IF
IF(ALLOCATED(key_table))THEN
  DEALLOCATE(key_table)
END IF
ALLOCATE(arr_int_2d(1:SIZE(arr_in,1),1:SIZE(arr_in,2)))
ALLOCATE(key_table(1: SIZE(arr_in,1)))
arr_int_2d(:,:) = arr_in(:,:)
DO n = 1, SIZE(arr_in,1)
  key_table(n) = n
END DO
END SUBROUTINE set_2d_int


SUBROUTINE set_1d_dp(arr_in)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: arr_in(:)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
choice = 3      ! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
IF(ALLOCATED(arr_dp_1d))THEN
  DEALLOCATE(arr_dp_1d)
END IF
IF(ALLOCATED(key_table))THEN
  DEALLOCATE(key_table)
END IF
ALLOCATE(arr_dp_1d(1:SIZE(arr_in,1)))
ALLOCATE(key_table(1:SIZE(arr_in,1)))
arr_dp_1d(:) = arr_in(:)
DO n = 1, SIZE(arr_in,1)
  key_table(n) = n
END DO
END SUBROUTINE set_1d_dp


SUBROUTINE set_2d_dp(arr_in)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: arr_in(:,:)
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
choice = 4      ! 1: int 1d, 2: int 2d, 3: dp 1d, 4: dp 2d
IF(ALLOCATED(arr_dp_2d))THEN
  DEALLOCATE(arr_dp_2d)
END IF
IF(ALLOCATED(key_table))THEN
  DEALLOCATE(key_table)
END IF
ALLOCATE(arr_dp_2d(1:SIZE(arr_in,1), 1:SIZE(arr_in,2)))
ALLOCATE(key_table(1:SIZE(arr_in,1)))
arr_dp_2d(:,:) = arr_in(:,:)
DO n = 1, SIZE(arr_in,1)
  key_table(n) = n
END DO
END SUBROUTINE set_2d_dp




