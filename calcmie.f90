subroutine calcmie(nWav,lambda,nrad,krad,a,Q)
  use constants_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,lambda(nWav),nrad(nWav),krad(nWav)
  double precision,intent(out) :: Q(nWav)
  double precision :: x,qext,qsca,gg
  double complex :: index
  integer :: i

  do i=1,nWav
     x = twopi*a/lambda(i)
     index = cmplx(nrad(i),krad(i))
     call bhmie(x,index,qext,qsca,gg)
     Q(i) = qext-qsca
     if (Q(i) .le. 0.0d0) Q(i) = 0.0d0
  end do

end subroutine calcmie
