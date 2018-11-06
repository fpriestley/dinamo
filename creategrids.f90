subroutine createHgrid(nEnthalpy,Tgrid,Tmid,Hgrid,Hmid,nTemp,H_T,T_H)

  implicit none

  integer,intent(in) ::nEnthalpy,nTemp
  double precision,intent(in) :: Tgrid(0:nEnthalpy),Tmid(nEnthalpy)
  double precision,intent(in) :: H_T(nTemp),T_H(nTemp)
  double precision,intent(out) :: Hgrid(0:nEnthalpy),Hmid(nEnthalpy)
  integer :: i

  do i=0,nEnthalpy
     call hinterp(Tgrid(i),Hgrid(i),T_H,H_T,nTemp)
  end do

  do i=1,nEnthalpy
     call hinterp(Tmid(i),Hmid(i),T_H,H_T,nTemp)
  end do

end subroutine createHgrid

subroutine createTgrid(nEnthalpy,Tgrid,Tmid,Tmin,Tmax)

  implicit none

  integer,intent(in) :: nEnthalpy
  double precision,intent(in) :: Tmin,Tmax
  double precision,intent(out) :: Tgrid(0:nEnthalpy),Tmid(nEnthalpy)
  double precision :: dT
  integer :: i

  dT = (Tmax-Tmin)/nEnthalpy

  do i=0,nEnthalpy
     Tgrid(i) = Tmin + i*dT
  end do

  do i=1,nEnthalpy
     Tmid(i) = 0.5d0*(Tgrid(i)+Tgrid(i-1))
  end do

end subroutine createTgrid

subroutine createsizegrid(nSizes,agrid,f_a,a_g,N_a)

  implicit none

  integer,intent(in) :: nSizes
  double precision,intent(in) :: agrid(0:nSizes),f_a(0:nSizes)
  double precision,intent(out) :: a_g(nSizes),N_a(nSizes)
  integer :: i

  do i=1,nSizes
     a_g(i) = 0.5d0*(agrid(i)+agrid(i-1))
     N_a(i) = 0.5d0*(f_a(i)+f_a(i-1))*(agrid(i)-agrid(i-1))
  end do

  N_a = 1d-4*N_a ! a in um, f in cm-1

end subroutine createsizegrid

subroutine createQ_a(nWav,a,Q_a,nQsize,Qsize,Qdata)

  implicit none

  integer,intent(in) :: nWav,nQsize
  double precision,intent(in) :: a,Qsize(nQsize),Qdata(nQsize,nWav)
  double precision,intent(out) :: Q_a(nWav)
  integer :: i

  do i=1,nWav
     call interpolate(a,Q_a(i),Qsize,Qdata(:,i),nQsize)
  end do

end subroutine createQ_a

subroutine createH_a(nTemp,H_T,H_a,a_g,rho_g)
  use constants_mod

  implicit none

  integer,intent(in) :: nTemp
  double precision,intent(in) :: a_g,rho_g,H_T(nTemp)
  double precision,intent(out) :: H_a(nTemp)
  integer :: i

  do i=1,nTemp
     H_a(i) = fourthirdpi*(1d-4*a_g)**3*rho_g*H_T(i)
  end do

end subroutine createH_a

subroutine createparticle(nEner,T,E,P)
  use constants_mod

  implicit none

  integer,intent(in) :: nEner
  double precision,intent(in) :: T
  double precision,intent(out) :: E(nEner),P(nEner)
  double precision :: Emax,dE,kT
  double precision :: boltzdist,integrate
  integer :: i

  kT = k_b*T

  Emax = 5.0d0*kT
  dE = Emax/nEner

  do i=1,nEner
     E(i) = i*dE
     P(i) = boltzdist(E(i),T)
  end do

end subroutine createparticle

subroutine mrndist(nSizes,amin,amax,n,agrid,f_a,N_g)

  implicit none

  integer,intent(in) :: nSizes
  double precision,intent(in) :: amin,amax,n,N_g
  double precision,intent(out) :: agrid(0:nSizes),f_a(0:nSizes)
  integer :: i
  double precision :: da,fda,integrate

  da = (amax/amin)**(1.0d0/nSizes)

  do i=0,nSizes
     agrid(i) = amin * da**i
     f_a(i) = agrid(i)**(-n)
  end do

  fda = integrate(f_a,agrid,nSizes+1,agrid(0),agrid(nSizes))
  f_a = 1d4*N_g*f_a/fda

end subroutine mrndist
