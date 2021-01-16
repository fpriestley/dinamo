subroutine readinput(runlabel,nSpecies,Nelec,N_H,Telec,T_atom,atomtype,radtype,radc1,radc2,T_cmb)

  implicit none

  character(len=50),intent(out) :: runlabel
  character(len=20),intent(out) :: radtype
  integer,intent(out) :: nSpecies
  double precision,intent(out) :: Nelec,N_H,Telec,T_atom,radc1,radc2,T_cmb
  character(len=2),intent(out) :: atomtype

  open(unit=1,file='input/input.dat',status='old')

  read(1,*) runlabel
  read(1,*) radtype,radc1,radc2
  read(1,*) T_cmb
  read(1,*) Nelec
  read(1,*) N_H
  read(1,*) Telec
  read(1,*) T_atom
  read(1,*) atomtype
  read(1,*) nSpecies

  close(unit=1)

end subroutine readinput

subroutine readgrainfiles(nSpecies,gfile)

  implicit none

  integer,intent(in) :: nSpecies
  character(len=100),intent(out) :: gfile(nSpecies)
  integer :: i

  open(unit=1,file='input/input.dat',status='old')

  do i=1,9
     read(1,*)
  end do

  do i=1,nSpecies
     read(1,*) gfile(i)
  end do

  close(unit=1)

end subroutine readgrainfiles

subroutine readprop(gfile,nkfile,rho,graintype,nSizes,amin,amax,nmrn,N_g,Tsub)

  implicit none

  character(len=100),intent(in) :: gfile
  character(len=100),intent(out) :: nkfile
  character(len=1),intent(out) :: graintype
  double precision,intent(out) :: rho,amin,amax,N_g,nmrn,Tsub
  integer,intent(out) :: nSizes

  open(unit=1,file='input/'//trim(gfile),status='old')

  read(1,*) graintype
  read(1,*) rho
  read(1,*) nkfile
  read(1,*) nSizes
  read(1,*) amin,amax
  read(1,*) nmrn
  read(1,*) N_g
  read(1,*) Tsub

  close(unit=1)

end subroutine readprop

subroutine readnk(nkfile,nWav,lambda,nrad,krad)

  implicit none

  character(len=100),intent(in) :: nkfile
  integer,intent(in) :: nWav
  double precision,intent(out) :: lambda(nWav),nrad(nWav),krad(nWav)
  character(len=3) :: type
  integer :: i

  type = nkfile(len(trim(nkfile))-2:len(trim(nkfile)))

  open(unit=1,file='input/'//trim(nkfile),status='old')

  if (type .eq. 'lnk') then
     do i=1,nWav
        read(1,*) lambda(i),nrad(i),krad(i)
     end do
  else if (type .eq. 'vnk') then
     do i=nWav,1,-1
        read(1,*) lambda(i),nrad(i),krad(i)
        lambda(i) = 1.0d4/lambda(i)
     end do
  else
     do i=1,nWav
        lambda(i) = 0.0d0
        nrad(i) = 0.0d0
        krad(i) = 0.0d0
     end do
  end if

  close(unit=1)

end subroutine readnk
