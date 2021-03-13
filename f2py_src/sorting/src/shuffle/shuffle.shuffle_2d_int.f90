SUBROUTINE shuffle_2d_int(arr_in, loops, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: arr_in(:,:)
INTEGER(KIND=StandardInteger), INTENT(IN) :: loops
INTEGER(KIND=StandardInteger), INTENT(OUT) :: arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) :: n
INTEGER(KIND=StandardInteger) :: a = 0
INTEGER(KIND=StandardInteger) :: b = 0
INTEGER(KIND=StandardInteger) :: t(1::SIZE(arr_in,2))
!#############################################################
arr(:,:) = arr_in(:,:)
DO n=1, loops * SIZE(arr_in,1)
  a = b
  DO WHILE(a .NE. b)
    CALL randint(1, SIZE(arr,1), a)
    CALL randint(1, SIZE(arr,1), b)
  END DO
  t(:) = arr(a,:)
  arr(a,:) = arr(b,:)
  arr(b,:) = t(:)
END DO
END SUBROUTINE shuffle_2d_int