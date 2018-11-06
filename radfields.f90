! ISRF function from Camps+2015
! J in erg cm-3 s-1 sr-1, lambda in um
double precision function JMat(lambda)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda
  double precision :: Jbb

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

! Fit to synchrotron spectrum - power law with exponential cutoff
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

! Fit to Crab Nebula synchrotron spectrum - power law with exponential cutoff
! Jcrab in erg cm-3 s-1 sr-1, lambda in um, d in pc
double precision function Jcrab(lambda,d)
  use constants_mod

  implicit none

  double precision,intent(in) :: lambda,d
  double precision :: nu,lambdacm,nulnu

  lambdacm = 1d-4*lambda
  nu = c_s/lambdacm ! frequency in Hz

  nulnu = 8.0d24 * nu**1.075 * exp(-(nu/2.4d5)**0.1) ! luminosity in erg s-1
  Jcrab = nulnu/fourpi/(d*pc)**2 ! flux at distance d pc in erg cm-2 s-1
  Jcrab = Jcrab/(pi*lambdacm) ! mean intensity in erg cm-3 s-1 sr-1

end function Jcrab
