SUBROUTINE sort_1d_dp(arr_in, ascending_in, arr)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) :: arr_in(:)
LOGICAL, OPTIONAL, INTENT(IN) :: ascending_in
REAL(kind=DoubleReal), INTENT(OUT) :: arr(1:SIZE(arr_in,1))
!#############################################################
REAL(kind=DoubleReal) :: temp
LOGICAL :: ascending = .TRUE.
INTEGER(KIND=StandardInteger) :: n, m
LOGICAL :: sorted = .FALSE.
!#############################################################
CALL CPU_TIME(start_time)

arr = arr_in

IF(PRESENT(ascending_in))THEN
  ascending = ascending_in
END IF

IF(ascending)THEN
  m = SIZE(arr,1)-1
  DO WHILE(.NOT. sorted)
    sorted = .TRUE.
    DO n = 1, m
      IF(arr(n) .GT. arr(n+1))THEN
        temp = arr(n)
        arr(n) = arr(n+1)
        arr(n+1) = temp
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
      IF(arr(n) .LT. arr(n+1))THEN
        temp = arr(n)
        arr(n) = arr(n+1)
        arr(n+1) = temp
        sorted = .FALSE.
      END IF
    END DO
    m = m - 1
  END DO    
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE sort_1d_dp