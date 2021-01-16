subroutine countlines(filename,n)

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(out) :: n
  integer :: io

  open(unit=99,file=filename,status='old')

  n = 0

  do
     read(99,*,iostat=io)
     if (io .ne. 0) exit
     n = n + 1
  end do

  close(unit=99)

end subroutine countlines
