SUBROUTINE sort_1d_dp_asc(arr_in, arr_out)
!#############################################################
REAL(kind=DoubleReal), INTENT(IN) ::         arr_in(:)
REAL(kind=DoubleReal), INTENT(OUT) ::        arr_out(1:SIZE(arr_in,1))
!#############################################################
INTEGER(KIND=StandardInteger) ::             n, m, st, kt
REAL(kind=DoubleReal) ::                     temp
LOGICAL ::                                   sort
!#############################################################
CALL make_keytable(SIZE(arr_in,1))
arr_out = arr_in
IF(SIZE(arr_in,1) .LE. mb_breakpoint)THEN
  m = SIZE(arr_out,1)-1
  sort = .TRUE.
  DO WHILE(sort)
    sort = .FALSE.
    DO n = 1, m
      IF(arr_out(n) .GT. arr_out(n+1))THEN
        temp = arr_out(n)
        arr_out(n) = arr_out(n+1)
        arr_out(n+1) = temp
        sort = .TRUE.
        ! Key table
        kt = keytable(n)
        keytable(n) = keytable(n+1)
        keytable(n+1) = kt        
      END IF
    END DO
    m = m - 1
  END DO
ELSE
  arr_out = arr_in
  CALL msort_1d_dp_asc(arr_out, keytable)
END IF
END SUBROUTINE sort_1d_dp_asc

RECURSIVE SUBROUTINE msort_1d_dp_asc(arr, ktarr)
!#############################################################
REAL(kind=DoubleReal), INTENT(INOUT) :: arr(:)
INTEGER(KIND=StandardInteger), INTENT(INOUT) :: ktarr(:)
!#############################################################
INTEGER(KIND=StandardInteger) :: lhs_len, rhs_len
INTEGER(KIND=StandardInteger) :: lhs_key, rhs_key
INTEGER(KIND=StandardInteger) :: m_key
REAL(kind=DoubleReal) :: merged(1:SIZE(arr, 1))
REAL(kind=DoubleReal) :: lhs(1:SIZE(arr,1) / 2)
REAL(kind=DoubleReal) :: rhs(1:SIZE(arr,1) - SIZE(arr,1) / 2)
INTEGER(KIND=StandardInteger) :: lhs_kt(1:SIZE(ktarr,1) / 2)
INTEGER(KIND=StandardInteger) :: rhs_kt(1:SIZE(ktarr,1) - SIZE(ktarr,1) / 2)
!#############################################################
IF(SIZE(arr,1) .GT. 1)THEN
  lhs_len = SIZE(arr,1) / 2
  rhs_len = SIZE(arr,1) - lhs_len
  
  CALL msort_1d_dp_asc(arr(1:lhs_len), ktarr(1:lhs_len))
  CALL msort_1d_dp_asc(arr(lhs_len+1:SIZE(arr,1)),ktarr(lhs_len+1:SIZE(arr,1)))
  
  lhs_key = 1
  rhs_key = 1
  m_key = 1
  
  lhs(:) = arr(1:lhs_len)
  rhs(:) = arr(lhs_len+1:SIZE(arr,1))  
  lhs_kt(:) = ktarr(1:lhs_len)
  rhs_kt(:) = ktarr(lhs_len+1:SIZE(ktarr,1))  
    
  IF(lhs(lhs_len) .LT. rhs(1))THEN
    arr(1:lhs_len) = lhs(:)
    arr(lhs_len+1:) = rhs(:) 
    ktarr(1:lhs_len) = lhs_kt(:)
    ktarr(lhs_len+1:) = rhs_kt(:)
  ELSE
    DO WHILE (lhs_key .LE. lhs_len .AND. rhs_key .LE. rhs_len)
      IF(lhs(lhs_key) .LE. rhs(rhs_key))THEN
        arr(m_key) = lhs(lhs_key)
        ktarr(m_key) = lhs_kt(lhs_key)
        lhs_key = lhs_key + 1
      ELSE
        arr(m_key) = rhs(rhs_key)
        ktarr(m_key) = rhs_kt(rhs_key)
        rhs_key = rhs_key + 1
      END IF
      m_key = m_key + 1
    END DO
    IF(lhs_key .GT. lhs_len)THEN
      arr(m_key:) = rhs(rhs_key:rhs_len)
      ktarr(m_key:) = rhs_kt(rhs_key:rhs_len)
    ELSE
      arr(m_key:) = lhs(lhs_key:lhs_len)
      ktarr(m_key:) = lhs_kt(lhs_key:lhs_len)
    END IF
  END IF  
END IF
END SUBROUTINE msort_1d_dp_asc