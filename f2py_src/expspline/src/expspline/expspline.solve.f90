SUBROUTINE normal_solve (a, b, x, flag)

IMPLICIT NONE

REAL(kind=DoubleReal), INTENT(IN) :: a(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: b(:)
REAL(kind=DoubleReal), INTENT(OUT) :: x(1:SIZE(a,1))
INTEGER(kind=StandardInteger), INTENT(OUT) :: flag


REAL(kind=DoubleReal) :: ata(SIZE(a,1),SIZE(a,1))
REAL(kind=DoubleReal) :: ata_c(SIZE(a,1),SIZE(a,1))
REAL(kind=DoubleReal) :: atb(SIZE(a,1))

flag = 0

ata(:,:) = matmul ( transpose ( a(:,:) ), a(:,:) )
atb(:) = matmul ( transpose ( a(:,:) ), b(:) )

call r8mat_cholesky_factor (ata, ata_c, flag )

if ( flag /= 0 ) then
  return
end if

call r8mat_cholesky_solve (ata_c, atb, x)

return
END SUBROUTINE normal_solve





subroutine r8mat_cholesky_factor (a, c, flag)

IMPLICIT NONE

REAL(kind=DoubleReal), INTENT(IN) :: a(:,:)
REAL(kind=DoubleReal), INTENT(OUT) :: c(1:SIZE(a,1),1:SIZE(a,1))
INTEGER(kind=StandardInteger), INTENT(OUT) :: flag

INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: i
INTEGER(kind=StandardInteger) :: j
REAL(kind=DoubleReal) :: sum2

  n = SIZE(a,1)
  flag = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end subroutine r8mat_cholesky_factor









subroutine r8mat_cholesky_solve (l, b, x )

IMPLICIT NONE

REAL(kind=DoubleReal), INTENT(IN) :: l(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: b(:)
REAL(kind=DoubleReal), INTENT(OUT) :: x(1:SIZE(l,1))
INTEGER(kind=StandardInteger) :: n
!
!  Solve L * y = b.
!
  call r8mat_l_solve (l, b, x )
!
!  Solve L' * x = y.
!
  call r8mat_lt_solve (l, x, x )

  return
end subroutine r8mat_cholesky_solve



subroutine r8mat_l_solve (a, b, x )

IMPLICIT NONE

REAL(kind=DoubleReal), INTENT(IN) :: a(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: b(:)
REAL(kind=DoubleReal), INTENT(OUT) :: x(1:SIZE(a,1))

INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: i 

!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end


subroutine r8mat_lt_solve (a, b, x )

IMPLICIT NONE

REAL(kind=DoubleReal), INTENT(IN) :: a(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: b(:)
REAL(kind=DoubleReal), INTENT(OUT) :: x(1:SIZE(a,1))

INTEGER(kind=StandardInteger) :: n
INTEGER(kind=StandardInteger) :: i 
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end