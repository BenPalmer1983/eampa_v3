SUBROUTINE inv(A_in, A_inverse)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: A_in(:,:)
REAL(kind=DoubleReal), INTENT(OUT) :: A_inverse(1:SIZE(A_in,1),1:SIZE(A_in,2))
!############################################################
REAL(kind=DoubleReal) :: A(1:SIZE(A_in,1),1:SIZE(A_in,2))
REAL(kind=DoubleReal) :: L(1:SIZE(A_in,1),1:SIZE(A_in,2))
REAL(kind=DoubleReal) :: U(1:SIZE(A_in,1),1:SIZE(A_in,2))
REAL(kind=DoubleReal) :: b(1:SIZE(A_in,1))
REAL(kind=DoubleReal) :: d(1:SIZE(A_in,1))
REAL(kind=DoubleReal) :: x(1:SIZE(A_in,1))
REAL(kind=DoubleReal) :: coeff
INTEGER(kind=StandardInteger) :: i, j, k, n
!############################################################

! SIZE
n = SIZE(A_in,1)

! INITIALISE MATRICES
A = A_in
L = 0.0D0
U = 0.0D0
b = 0.0D0
A_inverse = 0.0D0

! MAKE LOWER MATRIX
DO k=1,n-1
  DO i = k+1,n
    coeff = A(i,k) /A(k,k)
    L(i,k) = coeff
    DO j=k+1,n
      A(i,j) = A(i,j) - coeff * A(k,j)
    END DO
  END DO
END DO

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
DO i=1,n
  L(i,i) = 1.0
END DO
! U matrix is the upper triangular part of A
DO j=1,n
  DO i=1,j
    U(i,j) = a(i,j)
  END DO
END DO

! Step 3: compute columns of the inverse matrix C
DO k=1,n
  b(k) = 1.0D0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  DO i=2,n
    d(i) = b(i)
    DO j = 1, i-1
      d(i) = d(i) - L(i,j) * d(j)
    END DO
  END DO
! Step 3b: Solve Ux=d using the back substitution
  x(n) = d(n) / U(n,n)
  DO i = n-1, 1,-1
    x(i) = d(i)
    DO j = n, i+1, -1
      x(i) = x(i) - U(i,j) * x(j)
    END DO
    x(i) = x(i)/u(i,i)
  END DO
! Step 3c: fill the solutions x(n) into column k of C
  DO i=1,n
    A_inverse(i,k) = x(i)
  END DO
  b(k)=0.0
END DO
END SUBROUTINE inv