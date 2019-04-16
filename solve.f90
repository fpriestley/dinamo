subroutine solveTdist(nEnthalpy,P,nTemp,H_a,T_H,Tmid,nWav,lambda,Q_a,Jrad,a_g,rho_g,grain)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nEnthalpy,nTemp,nWav
  character(len=1),intent(in) :: grain
  double precision,intent(in) :: a_g,rho_g,H_a(nTemp),T_H(nTemp),lambda(nWav),Q_a(nWav),Jrad(nWav)
  double precision,intent(out) :: P(nEnthalpy),Tmid(nEnthalpy)
  double precision,parameter :: minP = 1d-20
  double precision,parameter :: minTdiff = 1.0d-2
  double precision :: Tmax,Tmin,Tequib
  double precision :: Tgrid(0:nEnthalpy),Hgrid(0:nEnthalpy),Hmid(nEnthalpy)
  double precision :: A(nEnthalpy,nEnthalpy),B(nEnthalpy)
  double precision :: zetaE,zetaEH
  integer :: j,i,info

  Tmin = 0.0d0
  Tmax = T_H(nTemp)

  if (lgElec .or. lgHyd) then
     do i=1,nEner
        if (lgElec) EzetaE(i) = zetaE(Eelec(i),a_g,rho_g,grain)*Eelec(i)
        if (lgHyd) EzetaE_H(i) = zetaEH(E_H(i),a_g,atomtype)*E_H(i)
     end do
  end if

  call equibT(nWav,lambda,Q_a,Jrad,a_g,Tequib,rho_g)

  if (lgEquibT) then
     do i=1,nEnthalpy
        Tmid(i) = Tequib + 1.0d0*(i-1)
     end do
     P = 0.0d0
     P(1) = 1.0d0
  else
     do
        call createTgrid(nEnthalpy,Tgrid,Tmid,Tmin,Tmax)
        call createHgrid(nEnthalpy,Tgrid,Tmid,Hgrid,Hmid,nTemp,H_a,T_H)
        call creatematrix(nEnthalpy,A,Hgrid,Hmid,Tmid,nWav,lambda,Q_a,Jrad,a_g)
        A(nEnthalpy,:) = 1.0d0
        B = 0.0d0
        B(nEnthalpy) = 1.0d0
        call linsolve(nEnthalpy,A,B,P,info)
        if (Tgrid(nEnthalpy) .ge. T_H(nTemp)) P(nEnthalpy) = 0.
        if (info .gt. 0) then
           Tmax = 0.9*Tmax
        else
           do j=nEnthalpy,1,-1
              if (P(j) .ge. minP*maxval(P)) then
                 Tmax = Tgrid(j)
                 exit
              end if
           end do
           if (Tmax .lt. 2.0d0*Tequib) then
              Tmax = 2.0d0*Tequib
           end if
        end if
        if (abs(Tmax-Tgrid(nEnthalpy)) .lt. minTdiff*Tmax) exit
     end do
  end if

  do i=1,nEnthalpy
     if (P(i) .lt. 1.0d-99) P(i) = 0.0d0
  end do

end subroutine solveTdist

subroutine equibT(nWav,lambda,Q,J,a,T0,rho)
  use particle_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,rho,lambda(nWav),Q(nWav),J(nWav)
  double precision,intent(out) :: T0
  double precision,parameter :: minTdiff = 0.1d0
  double precision :: Tmax,Tmin,Told,netheat
  double precision :: coolingrate,heatingrate

  T0 = 1.0d1
  Tmax = 1.0d2
  Tmin = 0.1d0

  do
     netheat = heatingrate(nWav,lambda,Q,J,a,rho) - coolingrate(T0,nWav,lambda,Q,a)
     if (netheat .gt. 0.0d0) then
        Told = T0
        Tmin = Told
        T0 = (Told + Tmax)/2.0d0
     else if (netheat .lt. 0.0d0) then
        Told = T0
        Tmax = Told
        T0 = (Told + Tmin)/2.0d0
     end if
     if (abs(T0-Told) .lt. minTdiff) exit
  end do

end subroutine equibT
