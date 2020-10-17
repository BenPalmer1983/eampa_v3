SUBROUTINE fit(x, y, order, coefficients)
! Fits a polynomial of order to the points input
! Uses Vandermonde Matrix
    Implicit None  !Force declaration of all variables
! Declare variables
! In
Real(kind=DoubleReal), INTENT(IN) ::                   x(:)
Real(kind=DoubleReal), INTENT(IN) ::                   y(:)
Integer(kind=StandardInteger), INTENT(IN) ::           order
! Out
Real(kind=DoubleReal), INTENT(OUT) ::                  coefficients(1:(order+1))
! Private
Integer(kind=StandardInteger) ::                       k,col,row,exponentValue
Real(kind=DoubleReal) ::                               xMatrix(1:(order+1),1:(order+1))
Real(kind=DoubleReal) ::                               yMatrix(1:(order+1))
! Build Least Squares Fitting Vandermonde matrix
Do row=1,(order+1)
  Do col=1,(order+1)
    exponentValue = row+col-2
    xMatrix(row,col) = 0.0D0
    Do k=1,size(x,1)
      xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*x(k)**exponentValue
    End Do
  End Do
End Do
Do row=1,(order+1)
  exponentValue = row-1
  yMatrix(row) = 0.0D0
  Do k=1,size(x,1)
    yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*y(k)*x(k)**exponentValue
  End Do
End Do

CALL solve(xMatrix, yMatrix, coefficients)

END SUBROUTINE fit