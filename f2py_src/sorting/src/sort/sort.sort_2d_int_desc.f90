SUBROUTINE sort_2d_int_desc(arr_in, col, arr_out)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::         arr_in(:, :)
INTEGER(KIND=StandardInteger), INTENT(IN) :: col
INTEGER(KIND=StandardInteger), INTENT(OUT) ::        arr_out(1:SIZE(arr_in,1), 1:SIZE(arr_in,2))
!#############################################################
INTEGER(KIND=StandardInteger) ::                     n, m, st
INTEGER(KIND=StandardInteger) ::                     temp(1:SIZE(arr_in,2))
LOGICAL ::                                           sort
!#############################################################
arr_out = arr_in
IF(SIZE(arr_in,1) .LE. 100)THEN
  m = SIZE(arr_out,1)-1
  sort = .TRUE.
  DO WHILE(sort)
    sort = .FALSE.
    DO n = 1, m
      IF(arr_out(n, col) .LT. arr_out(n+1, col))THEN
        temp(:) = arr_out(n,:)
        arr_out(n,:) = arr_out(n+1,:)
        arr_out(n+1,:) = temp(:)
        sort = .TRUE.
      END IF
    END DO
    m = m - 1
  END DO
ELSE
  arr_out = arr_in
  CALL msort_2d_int_desc(arr_out)
END IF
END SUBROUTINE sort_2d_int_desc

RECURSIVE SUBROUTINE msort_2d_int_desc(arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: arr(:, :)
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
INTEGER(KIND=StandardInteger) :: lhs_key, rhs_key
INTEGER(KIND=StandardInteger) :: m_key
INTEGER(KIND=StandardInteger) :: merged(1:SIZE(arr, 1), 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: lhs(1:SIZE(arr,1) / 2, 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: rhs(1:SIZE(arr,1) - SIZE(arr,1) / 2, 1:SIZE(arr, 2))
!#############################################################
IF(SIZE(arr,1) .GT. 1)THEN
  lhs_len = SIZE(arr,1) / 2
  rhs_len = SIZE(arr,1) - lhs_len
  
  CALL msort_2d_int_desc(arr(1:lhs_len,:))
  CALL msort_2d_int_desc(arr(lhs_len+1:SIZE(arr,1),:))
  
  lhs_key = 1
  rhs_key = 1
  m_key = 1
  
  lhs(:,:) = arr(1:lhs_len,:)
  rhs(:,:) = arr(lhs_len+1:SIZE(arr,1),:)  
  
  IF(lhs(lhs_len, col) .GT. rhs(1, col))THEN
    arr(1:lhs_len,:) = lhs(:,:)
    arr(lhs_len+1:,:) = rhs(:,:)  
  ELSE
    DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
      IF(lhs(lhs_key, col) .GE. rhs(rhs_key, col))THEN
        arr(m_key,:) = lhs(lhs_key,:)
        lhs_key = lhs_key + 1
      ELSE
        arr(m_key,:) = rhs(rhs_key,:)
        rhs_key = rhs_key + 1
      END IF
      m_key = m_key + 1
    END DO
    IF(lhs_key .GT. lhs_len)THEN
      arr(m_key:,:) = rhs(rhs_key:rhs_len,:)
    ELSE
      arr(m_key:,:) = lhs(lhs_key:lhs_len,:)
    END IF
  END IF  
END IF
END SUBROUTINE msort_2d_int_desc