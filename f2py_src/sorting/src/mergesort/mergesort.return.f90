SUBROUTINE ms_return_1d_dp(arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(OUT) :: arr(:)
!#############################################################
arr(:) = arr_dp_1d(:)
END SUBROUTINE ms_return_1d_dp


SUBROUTINE ms_return_2d_dp(arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(OUT) :: arr(:,:)
!#############################################################
arr(:,:) = arr_dp_2d(:,:)
END SUBROUTINE ms_return_2d_dp


SUBROUTINE ms_return_1d_int(arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(OUT) :: arr(:)
!#############################################################
arr(:) = arr_int_1d(:)
END SUBROUTINE ms_return_1d_int


SUBROUTINE ms_return_2d_int(arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(OUT) :: arr(:,:)
!#############################################################
arr(:,:) = arr_int_2d(:,:)
END SUBROUTINE ms_return_2d_int






