subroutine linsolve(N,A,B,X,info)

  implicit none

  integer,intent(in) :: N
  integer,intent(out) :: info
  double precision,intent(in) :: A(N,N),B(N)
  double precision,intent(out) :: X(N)
  double precision :: tempA(N,N)
  integer :: ipiv(N)

  tempA = A
  X = B

  call dgesv(N,1,tempA,N,ipiv,X,N,info)

end subroutine linsolve
