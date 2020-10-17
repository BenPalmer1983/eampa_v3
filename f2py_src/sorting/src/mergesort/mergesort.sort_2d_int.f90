SUBROUTINE sort_2d_int()
!#############################################################
!#############################################################
CALL CPU_TIME(start_time)
IF(make_key_table)THEN
  IF(sort_ascending)THEN
    CALL msort_2d_int_asc_kt(arr_int_2d, key_table)
  ELSE
    CALL msort_2d_int_desc_kt(arr_int_2d, key_table)
  END IF
ELSE
  IF(sort_ascending)THEN
    CALL msort_2d_int_asc(arr_int_2d)
  ELSE
    CALL msort_2d_int_desc(arr_int_2d)
  END IF
END IF
CALL CPU_TIME(end_time)
time = end_time - start_time
END SUBROUTINE sort_2d_int



RECURSIVE SUBROUTINE msort_2d_int_asc(arr)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: arr(:,:)
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
  
  CALL msort_2d_int_asc(arr(1:lhs_len,:))
  CALL msort_2d_int_asc(arr(lhs_len+1:SIZE(arr,1),:))
  
  lhs_key = 1
  rhs_key = 1
  m_key = 1
  
  lhs(:,:) = arr(1:lhs_len,:)
  rhs(:,:) = arr(lhs_len+1:SIZE(arr,1),:)  
  
  IF(lhs(lhs_len, col) .LT. rhs(1, col))THEN
    arr(1:lhs_len,:) = lhs(:,:)
    arr(lhs_len+1:,:) = rhs(:,:)  
  ELSE
    DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
      IF(lhs(lhs_key, col) .LE. rhs(rhs_key, col))THEN
        arr(m_key, :) = lhs(lhs_key, :)
        lhs_key = lhs_key + 1
      ELSE
        arr(m_key, :) = rhs(rhs_key, :)
        rhs_key = rhs_key + 1
      END IF
      m_key = m_key + 1
    END DO
    IF(lhs_key .GT. lhs_len)THEN
      arr(m_key:, :) = rhs(rhs_key:rhs_len, :)
    ELSE
      arr(m_key:, :) = lhs(lhs_key:lhs_len, :)
    END IF
  END IF  
END IF
END SUBROUTINE msort_2d_int_asc




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


RECURSIVE SUBROUTINE msort_2d_int_asc_kt(arr, kt)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: arr(:,:)
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: kt(:)
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
INTEGER(KIND=StandardInteger) :: lhs_key, rhs_key
INTEGER(KIND=StandardInteger) :: m_key
INTEGER(KIND=StandardInteger) :: merged(1:SIZE(arr, 1), 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: lhs(1:SIZE(arr,1) / 2, 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: rhs(1:SIZE(arr,1) - SIZE(arr,1) / 2, 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: lhs_kt(1:SIZE(arr,1) / 2)
INTEGER(KIND=StandardInteger) :: rhs_kt(1:SIZE(arr,1) - SIZE(arr,1) / 2)
!#############################################################
IF(SIZE(arr,1) .GT. 1)THEN
  lhs_len = SIZE(arr,1) / 2
  rhs_len = SIZE(arr,1) - lhs_len
  
  CALL msort_2d_int_asc_kt(arr(1:lhs_len,:), kt(1:lhs_len))
  CALL msort_2d_int_asc_kt(arr(lhs_len+1:SIZE(arr,1),:), kt(lhs_len+1:SIZE(arr,1)))
  
  lhs_key = 1
  rhs_key = 1
  m_key = 1
  
  lhs(:,:) = arr(1:lhs_len,:)
  rhs(:,:) = arr(lhs_len+1:SIZE(arr,1),:)  
  lhs_kt(:) = kt(1:lhs_len)
  rhs_kt(:) = kt(lhs_len+1:SIZE(arr,1))  
  
  IF(lhs(lhs_len, col) .LT. rhs(1, col))THEN
    arr(1:lhs_len,:) = lhs(:,:)
    arr(lhs_len+1:,:) = rhs(:,:)  
    kt(1:lhs_len) = lhs_kt(:) 
    kt(lhs_len+1:) = rhs_kt(:)  
  ELSE
    DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
      IF(lhs(lhs_key, col) .LE. rhs(rhs_key, col))THEN
        arr(m_key,:) = lhs(lhs_key,:)
        kt(m_key) = lhs_kt(lhs_key) 
        lhs_key = lhs_key + 1
      ELSE
        arr(m_key,:) = rhs(rhs_key,:)
        kt(m_key) = rhs_kt(rhs_key) 
        rhs_key = rhs_key + 1
      END IF
      m_key = m_key + 1
    END DO
    IF(lhs_key .GT. lhs_len)THEN
      arr(m_key:,:) = rhs(rhs_key:rhs_len,:)
      kt(m_key:) = rhs_kt(rhs_key:rhs_len)
    ELSE
      arr(m_key:,:) = lhs(lhs_key:lhs_len,:)
      kt(m_key:) = lhs_kt(lhs_key:lhs_len)
    END IF
  END IF  
END IF
END SUBROUTINE msort_2d_int_asc_kt




RECURSIVE SUBROUTINE msort_2d_int_desc_kt(arr, kt)
!#############################################################
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: arr(:,:)
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: kt(:)
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
INTEGER(KIND=StandardInteger) :: lhs_key, rhs_key
INTEGER(KIND=StandardInteger) :: m_key
INTEGER(KIND=StandardInteger) :: merged(1:SIZE(arr, 1), 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: lhs(1:SIZE(arr,1) / 2, 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: rhs(1:SIZE(arr,1) - SIZE(arr,1) / 2, 1:SIZE(arr, 2))
INTEGER(KIND=StandardInteger) :: lhs_kt(1:SIZE(arr,1) / 2)
INTEGER(KIND=StandardInteger) :: rhs_kt(1:SIZE(arr,1) - SIZE(arr,1) / 2)
!#############################################################
IF(SIZE(arr,1) .GT. 1)THEN
  lhs_len = SIZE(arr,1) / 2
  rhs_len = SIZE(arr,1) - lhs_len
  
  CALL msort_2d_int_desc_kt(arr(1:lhs_len,:), kt(1:lhs_len))
  CALL msort_2d_int_desc_kt(arr(lhs_len+1:SIZE(arr,1),:), kt(lhs_len+1:SIZE(arr,1)))
  
  lhs_key = 1
  rhs_key = 1
  m_key = 1
  
  lhs(:,:) = arr(1:lhs_len,:)
  rhs(:,:) = arr(lhs_len+1:SIZE(arr,1),:)  
  lhs_kt(:) = kt(1:lhs_len)
  rhs_kt(:) = kt(lhs_len+1:SIZE(arr,1))  
  
  IF(lhs(lhs_len, col) .GT. rhs(1, col))THEN
    arr(1:lhs_len,:) = lhs(:,:)
    arr(lhs_len+1:,:) = rhs(:,:)  
    kt(1:lhs_len) = lhs_kt(:) 
    kt(lhs_len+1:) = rhs_kt(:)  
  ELSE
    DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
      IF(lhs(lhs_key, col) .GE. rhs(rhs_key, col))THEN
        arr(m_key,:) = lhs(lhs_key,:)
        kt(m_key) = lhs_kt(lhs_key) 
        lhs_key = lhs_key + 1
      ELSE
        arr(m_key,:) = rhs(rhs_key,:)
        kt(m_key) = rhs_kt(rhs_key) 
        rhs_key = rhs_key + 1
      END IF
      m_key = m_key + 1
    END DO
    IF(lhs_key .GT. lhs_len)THEN
      arr(m_key:,:) = rhs(rhs_key:rhs_len,:)
      kt(m_key:) = rhs_kt(rhs_key:rhs_len)
    ELSE
      arr(m_key:,:) = lhs(lhs_key:lhs_len,:)
      kt(m_key:) = lhs_kt(lhs_key:lhs_len)
    END IF
  END IF  
END IF
END SUBROUTINE msort_2d_int_desc_kt



