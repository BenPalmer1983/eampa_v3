SUBROUTINE apply_keytable_1d_dp(arr_in, arr_out)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::         arr_in(:)
REAL(kind=DoubleReal), INTENT(OUT) ::        arr_out(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
IF (SIZE(arr_in, 1) .EQ. SIZE(keytable, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr_out(n) = arr_in(keytable(n))
  END DO
ELSE  
  arr_out(:) = arr_in(:)
END IF
!#############################################################
END SUBROUTINE apply_keytable_1d_dp




SUBROUTINE apply_keytable_2d_dp(arr_in, arr_out)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: arr_in(:,:)
REAL(kind=DoubleReal), INTENT(OUT) :: arr_out(1:SIZE(arr_in,1), 1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) :: n
!#############################################################
IF (SIZE(keytable, 1) .EQ. SIZE(arr_in, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr_out(n, :) = arr_in(keytable(n), :)
  END DO
ELSE IF (SIZE(keytable, 1) .EQ. SIZE(arr_in, 2)) THEN
  DO n=1, SIZE(arr_in, 2)
    arr_out(:, n) = arr_in(:, keytable(n))
  END DO
ELSE  
  arr_out(:,:) = arr_in(:,:)
END IF
!#############################################################
END SUBROUTINE apply_keytable_2d_dp