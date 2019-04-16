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
