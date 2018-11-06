! particle range R (4/3 a rho) in g cm-2 as a function of E in eV
! Eq. 13 from Dwek & Smith 1996 (eq. 3 from Dwek 1987 if lgOldRe = true)
double precision function R_E(E,grain)
  use constants_mod
  use particle_mod

  implicit none

  double precision,intent(in) :: E
  character(len=1),intent(in) :: grain
  integer :: i
  double precision :: logE,logR

  if (lgOldRe) then
     logE = log10(E)
     logR = -8.15d0 + 0.502d0*logE + 0.1466d0*logE**2
     R_E = 10**logR
  else
     logE = log10(E)
     if (grain .eq. 'C') then
        logR = agra(0)
        do i=1,4
           logR = logR + agra(i)*logE**i
        end do
     else if (grain .eq. 'S') then
        logR = asil(0)
        do i=1,4
           logR = logR + asil(i)*logE**i
        end do
     end if
     R_E = 10**logR
  end if

end function R_E

! particle energy E in eV as a function of range in g cm-2
double precision function E_R(R,grain)
  use constants_mod
  use particle_mod

  implicit none

  double precision,intent(in) :: R
  character(len=1),intent(in) :: grain
  double precision,parameter :: minRdiff = 1.0d-2
  double precision :: logR,logE
  double precision :: a,b,c
  double precision :: Emin,Emax,Rmin,Rmax
  double precision :: E0,E1,R0,Rdiff
  double precision :: R_E

  if (lgOldRe) then
     logR = log10(R)
     a = 0.1466d0
     b = 0.502d0
     c = -8.150d0 - logR
     if (c .gt. 0.0d0) c = 0.0d0
     logE = 0.5d0*(-b+sqrt(b**2 - 4.0d0*a*c))/a
     E_R = 10**logE
  else
     Emin = 0.0d0
     Emax = 1.0d6
     Rmin = R_E(Emin,grain)
     Rmax = R_E(Emax,grain)
     E1 = 1.0d2
     do
        R0 = R_E(E1,grain)
        if (R0 .lt. R) then
           E0 = E1
           Emin = E0
           E1 = (E0 + Emax)/2.0d0
        else if (R0 .gt. R) then
           E0 = E1
           Emax = E0
           E1 = (E0 + Emin)/2.0d0
        end if
        if (E1 .eq. Emin) exit
        Rdiff = abs(R_E(E1,grain)-R)/R
        if (Rdiff .lt. minRdiff) exit
     end do
     E_R = E1
  end if

end function E_R

! fraction of energy transferred by an electron of energy E (erg) to a dust grain
! of radius a (um), density rho (g cm-3)
! From Dwek & Smith 1996 (Dwek 1987 if lgOldRe = true)
double precision function zetaE(E,a,rho,grain)
  use constants_mod

  implicit none

  double precision,intent(in) :: E,a,rho
  character(len=1),intent(in) :: grain
  double precision :: R_E,E_R
  double precision :: EeV,Estar,Rstar,Edash,Rdash
  double precision :: E1,E2,dE,R0,R1,R2

  if (lgOldRe) then
     EeV = E/eV ! convert to eV from erg
     Rstar = 4.0d0/3.0d0 * rho * 1d-4*a ! critical grain thickness in g cm-2
     Estar = E_R(Rstar,grain)
     if (EeV .le. Estar) then
        zetaE = 0.875d0
     else
        Rdash = R_E(EeV,grain) - Rstar
        Edash = E_R(Rdash,grain)/EeV
        if (Edash .lt. 0.125d0) then
           zetaE = 0.875d0
        else
           zetaE = 1.0d0 - Edash
        end if
     end if
  else
     E1 = E/eV ! convert to eV from erg
     R0 = 4.0d0/3.0d0 * rho * 1d-4*a ! grain thickness in g cm-2
     R1 = R_E(E1,grain) ! range of electron in grain
     if (R1 .le. R0) then
        zetaE = 1.0d0
     else
        R2 = R1 - R0 ! excess electron range
        E2 = E_R(R2,grain) ! leftover electron energy
        dE = E1 - E2
        zetaE = dE/E1
     end if
  end if

end function zetaE

! fraction of energy transferred by incident hydrogen nuclei
! E in erg, a in um
! From Dwek 1987
double precision function zetaEH(E,a,type)
  use constants_mod

  implicit none

  double precision,intent(in) :: E,a
  character(len=2),intent(in) :: type
  double precision :: Ekev,Estar

  Ekev = 1.0d-3*E/eV ! convert to keV from erg

  if (type .eq. 'H ') then
     Estar = 1.33d2*a
  else if (type .eq. 'He') then
     Estar = 2.22d2*a
  else if (type .eq. 'C ') then
     Estar = 6.55d2*a
  else if (type .eq. 'O ') then
     Estar = 6.55d2*a
  end if

  if (Ekev .le. Estar) then
     zetaEH = 1.0d0
  else
     zetaEH = Estar/Ekev
  end if

end function zetaEH

! boltzmann distribution for particle energies
! E in erg, T in K
double precision function boltzdist(E,T)
  use constants_mod

  implicit none

  double precision,intent(in) :: E,T
  double precision :: kT

  kT = k_b*T
  boltzdist = 2.0d0/sqrt(pi) / (kT)**1.5 * sqrt(E) * exp(-E/kT)

end function boltzdist

! velocity of electron with energy E
! E in erg, v in cm s-1
double precision function v_E(E)
  use constants_mod

  implicit none

  double precision,intent(in) :: E

  v_E = sqrt(2.0d0*E/m_e)

end function v_E

! velocity of atom/ion with energy E
! E in erg, v in cm s-1
double precision function v_p(E,type)
  use constants_mod

  implicit none

  double precision,intent(in) :: E
  character(len=2),intent(in) :: type
  double precision :: mass

  if (type .eq. 'H ') then
     mass = m_p
  else if (type .eq. 'He') then
     mass = 2.0d0*m_p
  else if (type .eq. 'C ') then
     mass = 12.0d0*m_p
  else if (type .eq. 'O ') then
     mass = 16.0d0*m_p
  end if

  v_p = sqrt(2.0d0*E/mass)

end function v_p
