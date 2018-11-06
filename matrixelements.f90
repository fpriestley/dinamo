subroutine creatematrix(nEnthalpy,A,Hgrid,H,T,nWav,lambda,Q,J,a_g)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nEnthalpy,nWav
  double precision,intent(in) :: Hgrid(0:nEnthalpy),H(nEnthalpy),T(nEnthalpy)
  double precision,intent(in) :: lambda(nWav),Q(nWav),J(nWav),a_g
  double precision,intent(out) :: A(nEnthalpy,nEnthalpy)
  integer :: f,i
  double precision :: Hwidth(nEnthalpy),dUdt
  double precision :: coolingrate,contheat,heating_fi,heating_Ni,contcool,cooling_fi,cooling_1i

  A = 0.0d0

  do i=1,nEnthalpy
     Hwidth(i) = Hgrid(i) - Hgrid(i-1)
  end do

  do f=1,nEnthalpy
     do i=1,nEnthalpy
        if (i .eq. f) then
           cycle
        else if (i .gt. f) then
           if (f .eq. 1) then
              if (lgContcool) then
                 A(f,i) = 0.0d0
              else
                 A(f,i) = cooling_1i(H(f),H(i),Hwidth(f),Hgrid(0),nWav,T(i),lambda,Q,a_g)
              end if
           else
              if (lgContcool) then
                 A(f,i) = 0.0d0
              else
                 A(f,i) = cooling_fi(H(f),H(i),Hwidth(f),nWav,T(i),lambda,Q,a_g)
              end if
           end if
           if (i .eq. f+1) then
              if (lgContcool) then
                 dUdt = contheat(nWav,lambda,Q,J,a_g,Hgrid(i),H(i)) - coolingrate(T(i),nWav,lambda,Q,a_g)
              else
                 dUdt = contheat(nWav,lambda,Q,J,a_g,Hgrid(i),H(i)) - contcool(nWav,lambda,Q,a_g,H(i),Hgrid(i-1),T(i))
              end if
              if (dUdt .lt. 0.0d0) then
                 A(f,i) = A(f,i) - dUdt/Hwidth(i)
              end if
           end if
        else if (i .lt. f) then
           if (f .eq. nEnthalpy) then
              A(f,i) = heating_Ni(H(f),H(i),Hwidth(f),Hgrid(f),nWav,J,lambda,Q,a_g)
           else
              A(f,i) = heating_fi(H(f),H(i),Hwidth(f),nWav,J,lambda,Q,a_g)
           end if
           if (f .eq. i+1) then
              if (lgContcool) then
                 dUdt = contheat(nWav,lambda,Q,J,a_g,Hgrid(i),H(i)) - coolingrate(T(i),nWav,lambda,Q,a_g)
              else
                 dUdt = contheat(nWav,lambda,Q,J,a_g,Hgrid(i),H(i)) - contcool(nWav,lambda,Q,a_g,H(i),Hgrid(i-1),T(i))
              end if
              if (dUdt .gt. 0.0d0) then
                 A(f,i) = A(f,i) + dUdt/Hwidth(i)
              end if
           end if
        end if
     end do
  end do

  do i=1,nEnthalpy
     do f=1,nEnthalpy
        if (f .eq. i) then
           cycle
        else
           A(i,i) = A(i,i) + A(f,i)
        end if
     end do
     A(i,i) = -A(i,i)
  end do

end subroutine creatematrix

! Continuous cooling rate
! lambda, a in um, cooling rate in erg s-1
double precision function coolingrate(T,nWav,lambda,Q,a)
  use constants_mod

  implicit none

  double precision,intent(in) :: T,a
  integer,intent(in) :: nWav
  double precision,intent(in) :: lambda(nWav),Q(nWav)
  double precision :: integrand(nWav)
  double precision :: Jbb,integrate
  integer :: i

  do i=1,nWav
     integrand(i) = Q(i)*1d-4*Jbb(lambda(i),T) ! convert Jbb to account for lambda in um
  end do

  coolingrate = fourpisq*(a*1d-4)**2 * integrate(integrand,lambda,nWav,lambda(1),lambda(nWav))

end function coolingrate

! Heating rate used to determine equilibrium temperature
! lambda, a in um, heating rate in erg s-1
double precision function heatingrate(nWav,lambda,Q,J,a)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,lambda(nWav),Q(nWav),J(nWav)
  double precision :: integrand(nWav),eintegrand(nEner)
  double precision :: integrate,v_E,v_p
  integer :: i

  do i=1,nWav
     integrand(i) = Q(i)*1d-4*J(i) ! convert J to account for lambda in um
  end do

  heatingrate = fourpisq*(a*1d-4)**2 * integrate(integrand,lambda,nWav,lambda(1),lambda(nWav))

  ! electron collisional heating contribution
  if (lgElec) then
     do i=1,nEner
        eintegrand(i) = EzetaE(i)*Pelec(i)*v_E(Eelec(i))
     end do
     heatingrate = heatingrate + pi*(a*1d-4)**2 * Nelec * integrate(eintegrand,Eelec,nEner,Eelec(1),Eelec(nEner))
  end if

  ! hydrogen collisional heating contribution
  if (lgHyd) then
     do i=1,nEner
        eintegrand(i) = EzetaE_H(i)*P_H(i)*v_p(E_H(i),atomtype)
     end do
     heatingrate = heatingrate + pi*(a*1d-4)**2 * N_H * integrate(eintegrand,E_H,nEner,E_H(1),E_H(nEner))
  end if

end function heatingrate

! Heating by low energy photons within temperature bin
! lambda, a in um, heating rate in erg s-1
double precision function contheat(nWav,lambda,Q,J,a,H_u,H_l)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,lambda(nWav),Q(nWav),J(nWav),H_u,H_l
  double precision :: integrand(nWav),dH,lambdamin,eintegrand(nEner),dE
  double precision :: integrate,v_E,v_p
  integer :: i

  do i=1,nWav
     integrand(i) = Q(i)*1d-4*J(i) ! convert J to account for lambda in um
  end do

  dH = H_u - H_l

  lambdamin = 1d4*hc/dH ! convert to um

  contheat = fourpisq*(a*1d-4)**2 * integrate(integrand,lambda,nWav,lambdamin,lambda(nWav))

  ! electron collisional heating contribution
  if (lgElec) then
     do i=1,nEner
        dE = EzetaE(i) ! energy in erg transferred by electron with energy E
        if (dE .lt. dH) then
           eintegrand(i) = dE*Pelec(i)*v_E(Eelec(i))
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     contheat = contheat + pi*(a*1d-4)**2 * Nelec * integrate(eintegrand,Eelec,nEner,Eelec(1),Eelec(nEner))
  end if

  ! hydrogen collisional heating contribution
  if (lgHyd) then
     do i=1,nEner
        dE = EzetaE_H(i) ! energy in erg transferred by proton with energy E
        if (dE .lt. dH) then
           eintegrand(i) = dE*P_H(i)*v_p(E_H(i),atomtype)
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     contheat = contheat + pi*(a*1d-4)**2 * N_H * integrate(eintegrand,E_H,nEner,E_H(1),E_H(nEner))
  end if

end function contheat

! Continuous cooling within temperature bin
! lambda, a in um, cooling rate in erg s-1
double precision function contcool(nWav,lambda,Q,a,H_u,H_l,T)
  use constants_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: a,lambda(nWav),Q(nWav),H_u,H_l,T
  double precision :: integrand(nWav),dH,lambdamin
  double precision :: integrate,Jbb
  integer :: i

  do i=1,nWav
     integrand(i) = Q(i)*1d-4*Jbb(lambda(i),T) ! convert J to account for lambda in um
  end do

  dH = H_u - H_l

  lambdamin = 1d4*hc/dH ! convert to um

  contcool = fourpisq*(a*1d-4)**2 * integrate(integrand,lambda,nWav,lambdamin,lambda(nWav))

end function contcool

! Discrete heating from i to f
! lambda, a in um, heating rate in s-1
double precision function heating_fi(H_f,H_i,H_fwidth,nWav,J,lambda,Q,a)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: H_f,H_i,H_fwidth,a
  double precision,intent(in) :: J(nWav),lambda(nWav),Q(nWav)
  double precision :: invHdiff,lambdatrans,Jtrans,Qtrans,dH
  double precision :: dE,eintegrand(nEner)
  double precision :: v_E,v_p,integrate
  integer :: i

  dH = H_f - H_i
  invHdiff = 1.0d0/dH
  lambdatrans = 1d4*hc*invHdiff ! lambda in um
  
  if ((lambdatrans .gt. lambda(nWav)) .or. (lambdatrans .lt. lambda(1))) then
     heating_fi = 0.0d0
  else
     call interpolate(lambdatrans,Jtrans,lambda,J,nWav)
     call interpolate(lambdatrans,Qtrans,lambda,Q,nWav)
     heating_fi = fourpisq*(a*1d-4)**2*Qtrans*Jtrans*hc*H_fwidth*invHdiff**3
  end if

  ! electron collisional heating contribution
  if (lgElec) then
     do i=1,nEner
        dE = EzetaE(i) ! energy in erg transferred by electron with energy E
        if (abs(dE-dH) .lt. 0.5d0*H_fwidth) then
           eintegrand(i) = Pelec(i)*v_E(Eelec(i))
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     heating_fi = heating_fi + pi*(a*1d-4)**2 * Nelec * integrate(eintegrand,Eelec,nEner,Eelec(1),Eelec(nEner))
  end if

  ! hydrogen collisional heating contribution
  if (lgHyd) then
     do i=1,nEner
        dE = EzetaE_H(i) ! energy in erg transferred by proton with energy E
        if (abs(dE-dH) .lt. 0.5d0*H_fwidth) then
           eintegrand(i) = P_H(i)*v_p(E_H(i),atomtype)
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     heating_fi = heating_fi + pi*(a*1d-4)**2 * N_H * integrate(eintegrand,E_H,nEner,E_H(1),E_H(nEner))
  end if

end function heating_fi

! Heating to temperature of highest enthalpy bin or higher
! lambda, a in um, heating rate in s-1
double precision function heating_Ni(H_f,H_i,H_fwidth,Hmax,nWav,J,lambda,Q,a)
  use constants_mod
  use particle_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: H_f,H_i,H_fwidth,Hmax,a
  double precision,intent(in) :: J(nWav),lambda(nWav),Q(nWav)
  double precision :: Hmin,lambda_i,rategtN,integrand(nWav),lambdacm(nWav)
  double precision :: heating_fi,integrate,v_E,v_p
  double precision :: dE,eintegrand(nEner)
  integer :: i

  Hmin = Hmax - H_i ! minimum energy for transitions beyond highest bin
  lambda_i = hc/Hmin ! in cm
  lambdacm = 1d-4*lambda ! convert to cm for integration

  do i=1,nWav
     integrand(i) = lambdacm(i)*Q(i)*J(i)
  end do

  rategtN = fourpisq*(a*1d-4)**2/hc * integrate(integrand,lambdacm,nWav,lambdacm(1),lambda_i)

  ! electron collisional heating contribution
  if (lgElec) then
     do i=1,nEner
        dE = EzetaE(i) ! energy in erg transferred by electron with energy E
        if (dE .ge. Hmin) then
           eintegrand(i) = Pelec(i)*v_E(Eelec(i))
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     rategtN = rategtN + pi*(a*1d-4)**2 * Nelec * integrate(eintegrand,Eelec,nEner,Eelec(1),Eelec(nEner))
  end if

  ! hydrogen collisional heating contribution
  if (lgHyd) then
     do i=1,nEner
        dE = EzetaE_H(i) ! energy in erg transferred by proton with energy E
        if (dE .ge. Hmin) then
           eintegrand(i) = P_H(i)*v_p(E_H(i),atomtype)
        else
           eintegrand(i) = 0.0d0
        end if
     end do
     rategtN = rategtN + pi*(a*1d-4)**2 * N_H * integrate(eintegrand,E_H,nEner,E_H(1),E_H(nEner))
  end if

  heating_Ni = rategtN + heating_fi(H_f,H_i,H_fwidth,nWav,J,lambda,Q,a)

end function heating_Ni

! Discrete cooling from i to f
! lambda, a in um, cooling rate in s-1
double precision function cooling_fi(H_f,H_i,H_fwidth,nWav,T,lambda,Q,a)
  use constants_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: H_f,H_i,H_fwidth,a,T
  double precision,intent(in) :: lambda(nWav),Q(nWav)
  double precision :: invHdiff,lambdatrans,Qtrans,dH
  double precision :: Jbb
  integer :: i

  dH = H_i - H_f
  invHdiff = 1.0d0/dH
  lambdatrans = 1d4*hc*invHdiff ! lambda in um
  
  if ((lambdatrans .gt. lambda(nWav)) .or. (lambdatrans .lt. lambda(1))) then
     cooling_fi = 0.0d0
  else
     call interpolate(lambdatrans,Qtrans,lambda,Q,nWav)
     cooling_fi = fourpisq*(a*1d-4)**2*Qtrans*Jbb(lambdatrans,T)*hc*H_fwidth*invHdiff**3
  end if

end function cooling_fi

! Cooling to temperature of lowest enthalpy bin or lower
! lambda, a in um, cooling rate in s-1
double precision function cooling_1i(H_f,H_i,H_fwidth,Hmin,nWav,T,lambda,Q,a)
  use constants_mod

  implicit none

  integer,intent(in) :: nWav
  double precision,intent(in) :: H_f,H_i,H_fwidth,Hmin,a,T
  double precision,intent(in) :: lambda(nWav),Q(nWav)
  double precision :: Hmax,lambda_i,ratelt1,integrand(nWav),lambdacm(nWav)
  double precision :: cooling_fi,integrate,Jbb
  integer :: i

  Hmax = H_i - Hmin ! maximum energy for transitions beyong lowest bin
  lambda_i = hc/Hmax ! in cm
  lambdacm = 1d-4*lambda ! convert to cm for integration

  do i=1,nWav
     integrand(i) = lambdacm(i)*Q(i)*Jbb(lambda(i),T)
  end do

  ratelt1 = fourpisq*(a*1d-4)**2/hc * integrate(integrand,lambdacm,nWav,lambdacm(1),lambda_i)

  cooling_1i = ratelt1 + cooling_fi(H_f,H_i,H_fwidth,nWav,T,lambda,Q,a)

end function cooling_1i
