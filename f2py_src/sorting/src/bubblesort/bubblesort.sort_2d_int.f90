SUBROUTINE sort_2d_int(arr_in, col_in, ascending_in, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: arr_in(:,:)
INTEGER(KIND=StandardInteger), OPTIONAL, INTENT(IN) :: col_in
LOGICAL, OPTIONAL, INTENT(IN) :: ascending_in
INTEGER(KIND=StandardInteger), INTENT(OUT) :: arr(1:SIZE(arr_in,1),1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) :: col = 1
LOGICAL :: ascending = .TRUE.
INTEGER(KIND=StandardInteger) :: temp(1:SIZE(arr_in,2))
INTEGER(KIND=StandardInteger) :: n, m
LOGICAL :: sorted = .FALSE.
!#############################################################
CALL CPU_TIME(start_time)

arr = arr_in

IF(PRESENT(col_in))THEN
  col = col_in
END IF

IF(PRESENT(ascending_in))THEN
  ascending = ascending_in
END IF

IF(ascending)THEN
  m = SIZE(arr,1)-1
  DO WHILE(.NOT. sorted)
    sorted = .TRUE.
    DO n = 1, m
      IF(arr(n, col) .GT. arr(n+1, col))THEN
        temp(:) = arr(n, :)
        arr(n,:) = arr(n+1,:)
        arr(n+1,:) = temp(:)
        sorted = .FALSE.
      END IF
    END DO
    m = m - 1
  END DO
ELSE
  m = SIZE(arr,1)-1
  DO WHILE(.NOT. sorted)
    sorted = .TRUE.
    DO n = 1, m
      IF(arr(n, col) .LT. arr(n+1, col))THEN
        temp(:) = arr(n, :)
        arr(n,:) = arr(n+1,:)
        arr(n+1,:) = temp(:)
        sorted = .FALSE.
      END IF
    END DO
    m = m - 1
  END DO    
END IF

CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE sort_2d_int