subroutine readQparams(filename,nQsize,nWav)

  implicit none

  character(len=*),intent(in) :: filename
  character(len=5) :: junk
  integer,intent(out) :: nQsize,nWav
  integer :: i

  open(unit=1,file=filename)

  do i=1,7
     read(1,*)
  end do

  read(1,*) junk,nQsize
  read(1,*) junk,nWav

  close(unit=1)

end subroutine readQparams

subroutine readQ(filename,nQsize,nWav,Q,Qsize,lambda,rho_g)

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(in) :: nQsize,nWav
  double precision,intent(out) :: Q(nQsize,nWav),Qsize(nQsize),lambda(nWav),rho_g
  double precision :: junk
  integer :: i,j,nskip

  open(unit=1,file=filename)

  do i=1,9
     read(1,*)
  end do

  read(1,*) nskip,rho_g

  do i=1,nskip
     read(1,*)
  end do

  do i=1,nQsize
     read(1,*)
     read(1,*) Qsize(i)
     read(1,*)
     if (i .eq. 1) then
        do j=1,nWav
           read(1,*) junk,lambda(j),Q(i,j)
        end do
     else
        do j=1,nWav
           read(1,*) junk,junk,Q(i,j)
        end do
     end if
  end do

  close(unit=1)

end subroutine readQ

subroutine readnTemp(filename,nTemp)

  implicit none

  character(len=*),intent(in) :: filename
  character(len=24) :: junk
  integer,intent(out) :: nTemp
  double precision :: skip

  open(unit=1,file=filename)

  read(1,*)
  read(1,*)
  read(1,*) junk,junk,junk,skip,skip,nTemp

  close(unit=1)

end subroutine readnTemp

subroutine readHT(filename,nTemp,T_H,H_T)

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(in) :: nTemp
  double precision,intent(out) :: T_H(nTemp),H_T(nTemp)
  integer :: i

  open(unit=1,file=filename)

  do i=1,3
     read(1,*)
  end do

  do i=1,nTemp
     read(1,*) T_H(i),H_T(i)
  end do

  close(unit=1)

end subroutine readHT

subroutine readnSizes(filename,nSizes)

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(out) :: nSizes
  integer :: ierr

  open(unit=1,file=filename)

  nSizes = 0

  do
     read(1,*,iostat=ierr)
     if (ierr .ne. 0) exit
     nSizes = nSizes + 1
  end do

  nSizes = nSizes - 5

  close(unit=1)

end subroutine readnSizes

subroutine readSizes(filename,nSizes,a,f)

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(in) :: nSizes
  double precision,intent(out) :: a(0:nSizes),f(0:nSizes)
  integer :: i

  open(unit=1,file=filename)

  do i=1,4
     read(1,*)
  end do

  do i=0,nSizes
     read(1,*) a(i),f(i)
  end do

  close(unit=1)

end subroutine
