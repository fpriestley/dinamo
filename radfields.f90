! ISRF function from Camps+2015 with extinction from Cardelli+1989
! J in erg cm-3 s-1 sr-1, lambda in um, av in mag
double precision function JMat(lambda,av)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,av
  double precision :: Jbb,extinct

  JMat = 0.0d0

  if (lambda .lt. 0.0912) then
     JMat = 0.0d0
  else if (lambda .lt. 0.110) then
     JMat = 3.0690d4*lambda**3.4172
  else if (lambda .lt. 0.134) then
     JMat = 1.627d1
  else if (lambda .lt. 0.250) then
     JMat = 0.566d0*lambda**(-1.6678)
  else if (lambda .ge. 0.250) then
     JMat = 1d-14*Jbb(lambda,7.5d3) + 1d-13*Jbb(lambda,4.0d3) &
          & + 4d-13*Jbb(lambda,3.0d3)
  end if

  JMat = JMat * extinct(lambda,av)

end function JMat

! Blackbody spectrum
! J in erg cm-3 s-1 sr-1, lambda in um, T in K
double precision function Jbb(lambda,T)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,T
  double precision :: lambdacm,invlambda,boltzfac

  Jbb = 0.

  lambdacm = lambda*1d-4 ! convert to cm
  
  invlambda = 1d0/lambdacm
  boltzfac = 1d0/(k_b*T)

  Jbb = 2d0*hc*c_s*invlambda**5 / (exp(hc*invlambda*boltzfac) - 1d0)

end function Jbb

! Power law
! Jpl in erg cm-3 s-1 sr-1, lambda in um
double precision function Jpl(lambda,C,a)

  implicit none

  double precision,intent(in) :: lambda,C,a

  Jpl = C*lambda**a

end function Jpl

! Fit to synchrotron spectrum from Cas A - power law with exponential cutoff
! Jsync in erg cm-3 s-1 sr-1, lambda in um
double precision function Jsync(lambda,C,E0)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,C,E0
  double precision :: E,lambdacm

  lambdacm = 1d-4*lambda
  E = hc/lambdacm/eV ! energy in eV

  Jsync = C * E**0.2 * exp(-(E/E0)**0.4)
  Jsync = Jsync/(pi*lambdacm)

end function Jsync

! Double power law fit to Crab Nebula synchrotron spectrum
! Jcrab in erg cm-3 s-1 sr-1, lambda in um, d in pc
double precision function Jcrab(lambda,d,av)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,d,av
  double precision :: nu,lambdacm,nulnu,nu0,c,ex1,ex2,wavlim

  lambdacm = 1d-4*lambda
  nu = c_s/lambdacm ! frequency in Hz
  nu0 = 6.2467e13 ! break frequency
  c = 7.665e36 ! scaling factor
  ex1 = 0.630 ! first power law exponent
  ex2 = 0.0696 ! second power law exponent
  wavlim = 0.0912 ! lyman limit in um

  if (nu .le. nu0) then
     nulnu = c*(nu/nu0)**ex1
  else
     nulnu = c*(nu/nu0)**ex2
  end if
  Jcrab = nulnu/fourpi/(d*pc)**2 ! flux at distance d pc in erg cm-2 s-1
  Jcrab = Jcrab/(pi*lambdacm) ! mean intensity in erg cm-3 s-1 sr-1

  if ((av .le. 0) .and. (lambda .lt. wavlim)) Jcrab = 0.0

end function Jcrab

! Cardelli+1989 extinction law
! wav in um, av in mag, extinct the reduction factor
double precision function extinct(lambda,av)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,av
  double precision,parameter :: rv = 3.1
  double precision :: x,ax,bx,y,fa,fb,amag

  x = 1./lambda ! convert to wavenumber

  if (x .le. 1.1) then
     ax = 0.574*x**1.61
     bx = -0.527*x*1.61
  else if (x .le. 3.3) then
     y = x - 1.82
     ax = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
     bx = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
  else if (x .le. 5.9) then
     ax = 1.752 - 0.316*x - 0.104/((x-4.67)**2 + 0.341)
     bx = -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263)
  else if (x .le. 8.0) then
     y = x - 5.9
     fa = -0.04473*y**2 - 0.009779*y**3
     fb = 0.2130*y**2 + 0.1207*y**3
     ax = 1.752 - 0.316*x - 0.104/((x-4.67)**2 + 0.341) + fa
     bx = -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263) + fb
  else if (x .le. 10.0) then
     y = x - 8.0
     ax = -1.073 - 0.628*y + 0.137*y**2 - 0.070*y**3
     bx = 13.670 + 4.257*y - 0.420*y**2 + 0.374*y**3
  else
     ax = 100.
     bx = 100.
  end if

  amag = (ax + bx/rv) * av

  extinct = 10**(-0.4*amag)

end function extinct
