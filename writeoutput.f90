subroutine openoutput(nSpecies,label)

  implicit none

  integer,intent(in) :: nSpecies
  character(len=*),intent(in) :: label

  open(unit=10,file='output/'//trim(label)//'.Tdist.out',status='replace')
  open(unit=11,file='output/'//trim(label)//'.emis.out',status='replace')
  open(unit=12,file='output/'//trim(label)//'.param.out',status='replace')

  write(10,'(I1)') nSpecies
  write(11,'(I1)') nSpecies

end subroutine openoutput

subroutine closeoutput()

  implicit none

  close(unit=10)
  close(unit=11)
  close(unit=12)

end subroutine closeoutput

subroutine writeTdist(unitno,nEnthalpy,a_g,Tmid,P)

  implicit none

  integer,intent(in) :: unitno,nEnthalpy
  double precision,intent(in) :: a_g,Tmid(nEnthalpy),P(nEnthalpy)
  integer :: i

  write(unitno,'(ES10.3)') a_g
  write(unitno,'(1000F10.3)') (Tmid(i),i=1,nEnthalpy)
  write(unitno,'(1000ES10.3)') (P(i),i=1,nEnthalpy)

end subroutine writeTdist

subroutine writeemis(unitno,nWav,a_g,Jgrain)

  implicit none

  integer,intent(in) :: unitno,nWav
  double precision,intent(in) :: a_g,Jgrain(nWav)
  integer :: i

  write(unitno,'(ES10.3)') a_g
  write(unitno,'(10000ES10.3)') (Jgrain(i),i=1,nWav)

end subroutine writeemis

subroutine writeparams(unitno,radtype,radc1,radc2,Nelec,N_H,Telec,T_atom,atomtype,nSpecies,T_cmb)

  implicit none

  integer,intent(in) :: unitno,nSpecies
  character(len=20),intent(in) :: radtype
  double precision,intent(in) :: radc1,radc2,Nelec,N_H,Telec,T_atom,T_cmb
  character(len=2),intent(out) :: atomtype

  write(unitno,"(I2,2X,'No. grain species')") nSpecies
  write(unitno,"(A20,2X,2(ES10.3,2X),'Radiation field')") radtype,radc1,radc2
  write(unitno,"(ES10.3,2X,'CMB temperature (K)')") T_cmb
  write(unitno,"(ES10.3,2X,'Electron density (cm-3)')") Nelec
  write(unitno,"(ES10.3,2X,'Hydrogen density (cm-3)')") N_H
  write(unitno,"(ES10.3,2X,'Electron temperature (K)')") Telec
  write(unitno,"(ES10.3,2X,'Atom/ion temperature (K)')") T_atom
  write(unitno,"(A2,2X,'Atom/ion type')") atomtype

end subroutine writeparams

subroutine writegrainparams(unitno,nkfile,rho,graintype,nSizes,amin,amax,nmrn,N_g,Tsub)

  implicit none

  character(len=100),intent(in) :: nkfile
  character(len=1),intent(in) :: graintype
  double precision,intent(in) :: rho,amin,amax,N_g,nmrn,Tsub
  integer,intent(in) :: nSizes,unitno

  write(unitno,*)
  write(unitno,"(A100,2X,'Optical properties')") nkfile
  write(unitno,"(A1,2X,'Grain type')") graintype
  write(unitno,"(ES10.3,2X,'Density (g cm-3)')") rho
  write(unitno,"(ES10.3,2X,'Grain number density (cm-3)')") N_g
  write(unitno,"(I4,2X,'No. grain sizes')") nSizes
  write(unitno,"(F4.1,2X,2(ES10.3,2X),'Size distribution (n, amin (um), amax (um))')") nmrn,amin,amax
  write(unitno,"(ES10.3,2X,'Grain sublimation temperature (K)')") Tsub

end subroutine writegrainparams
