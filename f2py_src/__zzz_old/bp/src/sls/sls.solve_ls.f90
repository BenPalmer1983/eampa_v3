SUBROUTINE sls_solve_ls(lower, upper, y, x)
!############################################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: lower
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: upper
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: y
REAL(kind=DoubleReal), DIMENSION(:), INTENT(OUT) :: x
!############################################################
REAL(kind=DoubleReal), DIMENSION(1:size(lower,1)) :: b
INTEGER(kind=StandardInteger) :: m_size
INTEGER(kind=StandardInteger) :: row, col, rowT, colT
REAL(kind=DoubleReal) :: sum_term
!############################################################

m_size = SIZE(lower,1)

b = 0.0D0

!# Solve Lb = y
  Do row = 1, m_size
    sum_term = 0.0D0
    Do col = 1, row
      If (col .LT. row) Then
        sum_term = sum_term + lower(row, col) * b(col)
      Else If (col .EQ. row) Then
        b(row) = (y(row) - sum_term) / lower(row, col)
      End If
    End Do
  End Do

  !# Solve Ux=b
  Do rowT = 1, m_size
    row = m_size - (rowT-1)
    sum_term = 0.0D0
    Do colT =1, m_size
      col = m_size - (colT-1)
      If (col .GT. row) Then
        sum_term = sum_term + upper(row, col) * x(col)
      Else If (col .EQ. row) Then
        x(row) = (b(row) - sum_term) / upper(row, col)
      End If
    End Do
  End Do

END SUBROUTINE sls_solve_ls