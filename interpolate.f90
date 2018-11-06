subroutine interpolate(x,y,xtab,ytab,n)

  implicit none

  integer,intent(in) :: n
  double precision,intent(in) :: x,xtab(n),ytab(n)
  double precision,intent(out) :: y
  integer :: i
  double precision :: dx,grad

  if (x .lt. xtab(1)) then
     grad = (ytab(2)-ytab(1))/(xtab(2)-xtab(1))
     dx = xtab(1)-x
     y = ytab(1) - grad*dx
  else if (x .gt. xtab(n)) then
     grad = (ytab(n)-ytab(n-1))/(xtab(n)-xtab(n-1))
     dx = x-xtab(n)
     y = ytab(n) + grad*dx
  else
     do i=1,n-1
        if (xtab(i+1) .ge. x) then
           grad = (ytab(i+1)-ytab(i))/(xtab(i+1)-xtab(i))
           dx = x-xtab(i)
           y = ytab(i) + grad*dx
           exit
        end if
     end do
  end if

end subroutine interpolate

subroutine hinterp(x,y,xtab,ytab,n)

  implicit none

  integer,intent(in) :: n
  double precision,intent(in) :: x,xtab(n),ytab(n)
  double precision,intent(out) :: y

  if (x .lt. xtab(1)) then
     y = ytab(1)*(x/xtab(1))**4
  else if (x .gt. xtab(n)) then
     y = ytab(n) + (x - xtab(n))*(ytab(n)-ytab(n-1))/(xtab(n)-xtab(n-1))
  else
     call interpolate(x,y,xtab,ytab,n)
  end if

end subroutine hinterp
