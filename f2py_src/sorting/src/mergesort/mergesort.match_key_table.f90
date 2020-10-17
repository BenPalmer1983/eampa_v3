SUBROUTINE match_key_table_1d_int_kt(arr_in, kt, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::       arr_in(:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::       kt(:)
INTEGER(KIND=StandardInteger), INTENT(OUT) ::      arr(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(kt, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n) = arr_in(kt(n))
  END DO
ELSE  
  arr(:) = arr_in(:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_1d_int_kt

SUBROUTINE match_key_table_1d_int(arr_in, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::       arr_in(:)
INTEGER(KIND=StandardInteger), INTENT(OUT) ::      arr(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(key_table, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n) = arr_in(key_table(n))
  END DO
ELSE  
  arr(:) = arr_in(:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_1d_int





SUBROUTINE match_key_table_2d_int_kt(arr_in, kt, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::       arr_in(:,:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::       kt(:)
INTEGER(KIND=StandardInteger), INTENT(OUT) ::      arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(kt, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n,:) = arr_in(kt(n),:)
  END DO
ELSE  
  arr(:,:) = arr_in(:,:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_2d_int_kt

SUBROUTINE match_key_table_2d_int(arr_in, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::       arr_in(:,:)
INTEGER(KIND=StandardInteger), INTENT(OUT) ::      arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(key_table, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n,:) = arr_in(key_table(n),:)
  END DO
ELSE  
  arr(:,:) = arr_in(:,:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_2d_int





SUBROUTINE match_key_table_1d_dp_kt(arr_in, kt, arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::               arr_in(:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::       kt(:)
REAL(kind=DoubleReal), INTENT(OUT) ::              arr(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(kt, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n) = arr_in(kt(n))
  END DO
ELSE  
  arr(:) = arr_in(:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_1d_dp_kt

SUBROUTINE match_key_table_1d_dp(arr_in, arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::               arr_in(:)
REAL(kind=DoubleReal), INTENT(OUT) ::              arr(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(key_table, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n) = arr_in(key_table(n))
  END DO
ELSE  
  arr(:) = arr_in(:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_1d_dp





SUBROUTINE match_key_table_2d_dp_kt(arr_in, kt, arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::               arr_in(:,:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::       kt(:)
REAL(kind=DoubleReal), INTENT(OUT) ::              arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(kt, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n,:) = arr_in(kt(n),:)
  END DO
ELSE  
  arr(:,:) = arr_in(:,:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_2d_dp_kt

SUBROUTINE match_key_table_2d_dp(arr_in, arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::       arr_in(:,:)
REAL(kind=DoubleReal), INTENT(OUT) ::      arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) ::                   n
!#############################################################
CALL CPU_TIME(start_time)
IF (SIZE(arr_in, 1) .EQ. SIZE(key_table, 1)) THEN
  DO n=1, SIZE(arr_in, 1)
    arr(n,:) = arr_in(key_table(n),:)
  END DO
ELSE  
  arr(:,:) = arr_in(:,:)
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE match_key_table_2d_dp



