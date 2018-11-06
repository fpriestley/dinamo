double precision function integrate(ytab,xtab,n,a,b)

  implicit none

  integer,intent(in) :: n
  double precision,intent(in) :: ytab(n),xtab(n),a,b
  double precision :: dx,ans,grad,yint
  integer :: i

  ans = 0.0d0

  do i=1,n-1
     if (xtab(i+1) .lt. a) then
        cycle
     else if (xtab(i) .gt. b) then
        exit
     else if ((xtab(i) .lt. a) .and. (xtab(i+1) .gt. a)) then
        dx = xtab(i+1)-a
        grad = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
        yint = ytab(i) + grad*(a-xtab(i))
        ans = ans + 0.5d0*(ytab(i+1)+yint)*dx
     else if ((xtab(i) .lt. b) .and. (xtab(i+1) .gt. b)) then
        dx = b - xtab(i)
        grad = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
        yint = ytab(i) + grad*dx
        ans = ans + 0.5d0*(yint+ytab(i))*dx
     else if ((xtab(i) .ge. a) .and. (xtab(i+1) .le. b)) then
        dx = xtab(i+1)-xtab(i)
        ans = ans + 0.5d0*(ytab(i+1)+ytab(i))*dx
     end if
  end do

  integrate = ans

end function integrate
