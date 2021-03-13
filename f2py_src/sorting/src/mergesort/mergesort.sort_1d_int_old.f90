SUBROUTINE sort_1d_int(arr_in, ascending_in, arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) :: arr_in(:)
LOGICAL, OPTIONAL, INTENT(IN) :: ascending_in
INTEGER(KIND=StandardInteger), INTENT(OUT) :: arr(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) :: temp
LOGICAL :: ascending = .TRUE.
INTEGER(KIND=StandardInteger) :: n, m
LOGICAL :: sorted = .FALSE.
!#############################################################
CALL CPU_TIME(start_time)
arr = arr_in
CALL msort_1d_int(arr, ascending_in)
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE sort_1d_int



RECURSIVE SUBROUTINE msort_1d_int (arr, ascending_in)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: arr(:)
LOGICAL, OPTIONAL, INTENT(IN) :: ascending_in
!#############################################################
LOGICAL :: ascending = .TRUE.
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
!#############################################################
IF(SIZE(arr,1) .GT. 1)THEN
  lhs_len = SIZE(arr,1) / 2
  rhs_len = SIZE(arr,1) - lhs_len
  CALL msort_1d_int(arr(1:lhs_len), ascending)
  CALL msort_1d_int(arr(lhs_len+1:SIZE(arr,1)), ascending)
  arr(1:SIZE(arr,1)) = msort_1d_int_merge(arr(1:lhs_len), arr(lhs_len+1:SIZE(arr,1)), ascending)
END IF
END SUBROUTINE msort_1d_int




FUNCTION msort_1d_int_merge(lhs, rhs, ascending) RESULT (merged)
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs(:)
INTEGER(KIND=StandardInteger) :: rhs(:)
LOGICAL :: ascending
INTEGER(KIND=StandardInteger) :: merged(1:SIZE(lhs, 1)+SIZE(rhs, 1))
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
INTEGER(KIND=StandardInteger) :: lhs_key, rhs_key
INTEGER(KIND=StandardInteger) :: m_key
!#############################################################
lhs_len = SIZE(lhs,1)
rhs_len = SIZE(rhs,1)
lhs_key = 1
rhs_key = 1
m_key = 1

IF(lhs(lhs_len) .LT. rhs(1))THEN
  merged(1:lhs_len) = lhs(:)
  merged(lhs_len+1:) = rhs(:)
ELSE
  DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
    IF(lhs(lhs_key) .LE. rhs(rhs_key))THEN
      merged(m_key) = lhs(lhs_key)
      lhs_key = lhs_key + 1
    ELSE
      merged(m_key) = rhs(rhs_key)
      rhs_key = rhs_key + 1
    END IF
    m_key = m_key + 1
  END DO
  IF(lhs_key .GT. lhs_len)THEN
    merged(m_key:) = rhs(rhs_key:rhs_len)
  ELSE
    merged(m_key:) = lhs(lhs_key:lhs_len)
  END IF
END IF

END FUNCTION msort_1d_int_merge






