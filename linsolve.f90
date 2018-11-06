subroutine linsolve(N,A,B,X)

  implicit none

  integer,intent(in) :: N
  double precision,intent(in) :: A(N,N),B(N)
  double precision,intent(out) :: X(N)
  double precision :: tempA(N,N)
  integer :: ipiv(N),info

  tempA = A
  X = B

  call dgesv(N,1,tempA,N,ipiv,X,N,info)

end subroutine linsolve
